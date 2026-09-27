// forged_dec_shim.cpp — pure-C ABI over Forged-lab's audited native DEC, so the
// Julia thesis experiments run REAL discrete exterior calculus, not hand-rolled
// (previously faked) operators.
//
// The metric Hodge 1-Laplacian uses the FEEC / WHITNEY route, which is the
// principled one on arbitrary (non-well-centered / sheared / obtuse) meshes:
//
//   * linarm::SimplicialComplex  — the explicit simplicial complex over the
//     actual mesh: signed boundary maps B_k (∂∂=0, verified via dd_max), Betti.
//   * linarm::galerkin_hodge_star — the Whitney FEEC mass matrices
//     M_k[σ,τ] = ∫⟨W_σ,W_τ⟩dV (Arnold–Falk–Winther; Hiptmair edge elements).
//     ALWAYS SPD and mesh-robust — unlike the diagonal circumcentric ★, which
//     goes singular/negative on non-well-centered triangles. This is why the
//     old 0.25·|e| floor + ★1=0 guards are gone: with the Galerkin mass there
//     is nothing to patch.
//
// For a 2D mesh the 2-form (top) Whitney mass is the robust diagonal 1/area, so
// the full Hodge-1  Δ₁ = δ₂d₁ + d₀δ₁  is assembled cleanly from M0,M1 (FEEC) and
// M2 = diag(1/area). The shim hands Julia the signed boundaries + SPD mass
// matrices; Julia solves the symmetric-definite generalized eigenproblem.
//
// Also exposes the audited 3D circumcentric star (build_hodge3d, the HKV-2017
// positivity-proven case) and the SU(2) Wilson-loop holonomy.

#include "geom/manifold.h"
#include "geom/dec/hodge.h"
#include "linarm/dec.hpp"
#include "linarm/hodge_galerkin.h"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>

extern "C" {
void geom_quat_exp_f32(const float *v3, float *q4, size_t n);
void geom_su2_wilson_f32(const float *edges, int k, float *W, float *trace, size_t n);
// MagNet directed-graph Laplacian (arXiv:2102.11391): complex Hermitian L^q
// whose U(1) phase Θ^q = 2πq(A_uv − A_vu) encodes edge DIRECTION. Hermitian ⇒
// real spectrum (no symmetrize/SVD hack). L = D_s − A_s∘e^{iΘ^q}.
void geom_magnetic_laplacian_f32(const float *A, int n, float q, float *L_re, float *L_im);
}

// ── Native NEON integer NTT (generic NTT-friendly prime) — exact, for Falcon
// mod q=12289. Cyclic; the negacyclic (xⁿ+1) ψ-weighting is done by the caller. ─
extern "C" {
void geom_ntt_twiddles(uint32_t root, uint32_t p, size_t n, uint32_t *tw, uint32_t *tw_inv);
void geom_ntt_forward(uint32_t *a, const uint32_t *tw, uint32_t p, size_t n);
void geom_ntt_inverse(uint32_t *a, const uint32_t *tw_inv, uint32_t p, uint32_t n_inv, size_t n);
}
extern "C" void forged_ntt_twiddles(uint32_t root, uint32_t p, size_t n, uint32_t *tw, uint32_t *twi) {
    geom_ntt_twiddles(root, p, n, tw, twi);
}
extern "C" void forged_ntt_forward(uint32_t *a, const uint32_t *tw, uint32_t p, size_t n) {
    geom_ntt_forward(a, tw, p, n);
}
extern "C" void forged_ntt_inverse(uint32_t *a, const uint32_t *twi, uint32_t p, uint32_t ni, size_t n) {
    geom_ntt_inverse(a, twi, p, ni, n);
}

// ── Directed-graph magnetic Laplacian L^q (n×n, row-major re/im) ────────────
// A: directed adjacency (n×n row-major). q: chirality / θ-angle. Direction
// lives in A−Aᵀ; q=0 collapses to the real symmetric graph Laplacian.
extern "C" int forged_magnetic_laplacian(const float *A, int n, float q,
                                         float *L_re, float *L_im) {
    if (!A || n <= 0) return -1;
    geom_magnetic_laplacian_f32(A, n, q, L_re, L_im);
    return 0;
}

// Densify a linarm signed Boundary (col-sparse, rows×cols) into row-major out.
static void densify(const linarm::Boundary &B, double *out) {
    std::memset(out, 0, sizeof(double) * static_cast<size_t>(B.rows) * B.cols);
    for (int j = 0; j < B.cols; ++j)
        for (auto &rs : B.col[static_cast<size_t>(j)])
            out[static_cast<size_t>(rs.first) * B.cols + j] = static_cast<double>(rs.second);
}

// ── FEEC / Whitney metric Hodge-1 building blocks for an ANY-dimension mesh ──
// pts: point-major nV*dim doubles. simps: nS top simplices, each `sv` ascending
// vertex indices (sv = topdim+1: 3 triangles 2D, 4 tets 3D, 5 pentatopes 4D).
// Query phase (B1==nullptr): writes nE,nF via out_nE/out_nF and dd2 (∂₁∂₂ max
// entry, must be 0), returns 0. Fill phase (row-major):
//   B1 [nV*nE]  signed ∂₁ (rows=verts, cols=edges)
//   B2 [nE*nF]  signed ∂₂ (rows=edges, cols=faces)
//   M0 [nV*nV]  Whitney 0-form (P1) mass, SPD
//   M1 [nE*nE]  Whitney 1-form (edge)  mass, SPD
//   M2 [nF*nF]  Whitney 2-form (face)  mass, SPD  ← makes the metric Hodge-1
//               mesh-robust in 3D/4D too (no diagonal circumcentric ★, which
//               breaks on right/obtuse simplices). All masses from
//               linarm::galerkin_hodge_star; ∂∂=0 from the explicit complex.
extern "C" int forged_feec_hodge1(const double *pts, int nV, int dim,
                                  const int *simps, int nS, int sv,
                                  int *out_nE, int *out_nF, long long *out_dd2,
                                  double *B1, double *B2,
                                  double *M0, double *M1, double *M2) {
    if (!pts || !simps || nV <= 0 || nS <= 0 || sv < 3 || dim < 2) return -1;

    std::vector<std::vector<int>> simplices(static_cast<size_t>(nS));
    for (int s = 0; s < nS; ++s) {
        simplices[static_cast<size_t>(s)].assign(static_cast<size_t>(sv), 0);
        for (int j = 0; j < sv; ++j)
            simplices[static_cast<size_t>(s)][static_cast<size_t>(j)] = simps[s*sv + j];
    }
    linarm::SimplicialComplex C = linarm::SimplicialComplex::from_simplices(simplices);

    const int nE = C.count(1), nF = C.count(2);
    if (out_nE) *out_nE = nE;
    if (out_nF) *out_nF = nF;
    if (out_dd2) *out_dd2 = C.dd_max(2);     // ∂₁∂₂ max entry — correctness invariant (0)
    if (!B1) return 0;                       // query phase

    densify(C.boundary(1), B1);
    densify(C.boundary(2), B2);
    std::vector<double> coords(pts, pts + static_cast<size_t>(nV) * dim);
    std::vector<double> m0 = linarm::galerkin_hodge_star(C, coords, dim, 0);
    std::vector<double> m1 = linarm::galerkin_hodge_star(C, coords, dim, 1);
    std::vector<double> m2 = linarm::galerkin_hodge_star(C, coords, dim, 2);
    std::memcpy(M0, m0.data(), sizeof(double) * m0.size());
    std::memcpy(M1, m1.data(), sizeof(double) * m1.size());
    std::memcpy(M2, m2.data(), sizeof(double) * m2.size());
    return 0;
}

// Combinatorial Betti number β_k of the complex (nullity of the unit Hodge
// Laplacian) — the metric-independent ground truth for dim ker(Δ_k). Used to
// validate the FEEC harmonic count (Hodge theorem). simps: nS×sv vertex indices.
extern "C" int forged_betti(const int *simps, int nS, int sv, int k) {
    if (!simps || nS <= 0 || sv < 2) return -1;
    std::vector<std::vector<int>> simplices(static_cast<size_t>(nS));
    for (int s = 0; s < nS; ++s) {
        simplices[static_cast<size_t>(s)].assign(static_cast<size_t>(sv), 0);
        for (int j = 0; j < sv; ++j)
            simplices[static_cast<size_t>(s)][static_cast<size_t>(j)] = simps[s*sv + j];
    }
    return linarm::SimplicialComplex::from_simplices(simplices).betti(k);
}

using geom::Manifold;

// ── 3D circumcentric Hodge stars (audited build_hodge3d; HKV-2017 3D case) ──
extern "C" int forged_hodge3d(const float *pts, int nV,
                              const uint32_t *tets, int nTet,
                              int *out_nE, int *out_nF,
                              float *star0, float *star1, float *star2, float *star3,
                              uint32_t *edges, uint32_t *faces) {
    if (!pts || !tets || nV <= 0 || nTet <= 0) return -1;
    std::vector<float> coords(pts, pts + static_cast<size_t>(nV) * 3);
    Manifold m = Manifold::point_cloud_nd(std::move(coords), 3);
    std::vector<uint32_t> tflat(tets, tets + static_cast<size_t>(nTet) * 4);
    m.add_simplices(3, std::move(tflat));

    geom::Hodge3D H = geom::build_hodge3d(m);
    const int nE = static_cast<int>(H.edges.size());
    const int nF = static_cast<int>(H.faces.size());
    if (out_nE) *out_nE = nE;
    if (out_nF) *out_nF = nF;
    if (!star0) return 0;

    std::memcpy(star0, H.star0.data(), sizeof(float) * static_cast<size_t>(nV));
    std::memcpy(star1, H.star1.data(), sizeof(float) * static_cast<size_t>(nE));
    std::memcpy(star2, H.star2.data(), sizeof(float) * static_cast<size_t>(nF));
    std::memcpy(star3, H.star3.data(), sizeof(float) * static_cast<size_t>(nTet));
    for (int e = 0; e < nE; ++e) { edges[2*e] = H.edges[e][0]; edges[2*e+1] = H.edges[e][1]; }
    for (int f = 0; f < nF; ++f) {
        faces[3*f] = H.faces[f][0]; faces[3*f+1] = H.faces[f][1]; faces[3*f+2] = H.faces[f][2];
    }
    return 0;
}

// ── SU(2) Wilson-loop holonomy of a discrete connection ─────────────────────
// half_angle_vecs: k vectors of 3 floats; each is (θ_i/2)·n̂_i (geom_quat_exp →
// unit rotor). W4 ← ordered product (geom_su2_wilson_f32); trace ← tr(W).
extern "C" void forged_su2_wilson(const float *half_angle_vecs, int k,
                                  float *W4, float *trace) {
    std::vector<float> quats(static_cast<size_t>(k) * 4);
    for (int i = 0; i < k; ++i)
        geom_quat_exp_f32(half_angle_vecs + 3 * i, quats.data() + 4 * i, 1);
    geom_su2_wilson_f32(quats.data(), k, W4, trace, 1);
}
