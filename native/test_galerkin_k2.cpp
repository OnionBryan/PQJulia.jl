// test_galerkin_k2.cpp — guards the k=2 (Whitney 2-form / face) mass added to
// native linarm::galerkin_hodge_star. Two discriminating checks:
//   (1) 2D: the single 2-face of a triangle is the top form ⇒ M₂ = [1/area].
//   (2) 3D: a tetrahedron's 4-face mass M₂ is symmetric positive-definite
//       (Cholesky succeeds) — the property the FEEC Hodge-1 needs.
// Build: see build_shim.sh's compile flags (links libm5_math.a + the recompiled
// hodge_galerkin.o).  Returns 0 on success.

#include "linarm/dec.hpp"
#include "linarm/hodge_galerkin.h"

#include <cmath>
#include <cstdio>
#include <vector>

static int fails = 0;
static void check(bool c, const char *m) {
    std::printf("  %s  %s\n", c ? "ok  " : "FAIL", m);
    if (!c) ++fails;
}

// Cholesky succeeds ⇔ symmetric positive-definite (n×n row-major).
static bool is_spd(const std::vector<double> &A, int n) {
    std::vector<double> L((size_t)n * n, 0.0);
    for (int j = 0; j < n; ++j) {
        double s = A[(size_t)j * n + j];
        for (int k = 0; k < j; ++k) s -= L[(size_t)j * n + k] * L[(size_t)j * n + k];
        if (s <= 0.0) return false;
        L[(size_t)j * n + j] = std::sqrt(s);
        for (int i = j + 1; i < n; ++i) {
            double a = A[(size_t)i * n + j];
            for (int k = 0; k < j; ++k) a -= L[(size_t)i * n + k] * L[(size_t)j * n + k];
            L[(size_t)i * n + j] = a / L[(size_t)j * n + j];
        }
    }
    return true;
}

int main() {
    using linarm::SimplicialComplex;
    using linarm::galerkin_hodge_star;

    // (1) 2D triangle, area 1/2 ⇒ Whitney top-form mass = 1/area = 2.
    {
        SimplicialComplex C = SimplicialComplex::from_simplices({{0, 1, 2}});
        std::vector<double> coords = {0, 0, 1, 0, 0, 1};                 // 3 verts × dim 2
        auto M2 = galerkin_hodge_star(C, coords, 2, 2);
        check(C.count(2) == 1, "2D: one 2-face");
        check(M2.size() == 1 && std::fabs(M2[0] - 2.0) < 1e-9, "2D: M2 = 1/area = 2.0 (top form)");
    }

    // (2) 3D tetrahedron: 4-face Whitney 2-form mass is SPD.
    {
        SimplicialComplex C = SimplicialComplex::from_simplices({{0, 1, 2, 3}});
        std::vector<double> coords = {0, 0, 0,  1, 0.1, 0,  0.2, 1, 0.1,  0.1, 0.2, 1}; // 4×dim3
        auto M2 = galerkin_hodge_star(C, coords, 3, 2);
        const int nF = C.count(2);
        check(nF == 4, "3D: tet has 4 faces");
        // symmetry
        bool sym = true;
        for (int i = 0; i < nF; ++i)
            for (int j = 0; j < nF; ++j)
                if (std::fabs(M2[(size_t)i * nF + j] - M2[(size_t)j * nF + i]) > 1e-9) sym = false;
        check(sym, "3D: M2 symmetric");
        check(is_spd(M2, nF), "3D: M2 SPD (Cholesky succeeds)");
    }

    if (fails) { std::printf("\n%d k=2 GALERKIN CHECK(S) FAILED\n", fails); return 1; }
    std::printf("\nnative galerkin k=2 (Whitney 2-form mass): SOLID\n");
    return 0;
}
