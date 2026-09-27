#include <oqs/oqs.h>
#include <string.h>
/* Size accessors and thin wrappers so Julia can ccall liboqs without mirroring struct layouts. */
void *shim_sig_new(const char *n) { return OQS_SIG_new(n); }
size_t shim_sig_len(void *s, int which) { OQS_SIG *g = s; return which==0 ? g->length_public_key : which==1 ? g->length_secret_key : g->length_signature; }
int shim_sig_keypair(void *s, uint8_t *pk, uint8_t *sk) { return OQS_SIG_keypair(s, pk, sk); }
int shim_sig_sign(void *s, uint8_t *sig, size_t *siglen, const uint8_t *m, size_t mlen, const uint8_t *ctx, size_t ctxlen, const uint8_t *sk) {
  return ctxlen ? OQS_SIG_sign_with_ctx_str(s, sig, siglen, m, mlen, ctx, ctxlen, sk) : OQS_SIG_sign(s, sig, siglen, m, mlen, sk); }
int shim_sig_verify(void *s, const uint8_t *m, size_t mlen, const uint8_t *sig, size_t siglen, const uint8_t *ctx, size_t ctxlen, const uint8_t *pk) {
  return ctxlen ? OQS_SIG_verify_with_ctx_str(s, m, mlen, sig, siglen, ctx, ctxlen, pk) : OQS_SIG_verify(s, m, mlen, sig, siglen, pk); }
void *shim_kem_new(const char *n) { return OQS_KEM_new(n); }
size_t shim_kem_len(void *k, int which) { OQS_KEM *g = k; return which==0 ? g->length_public_key : which==1 ? g->length_secret_key : which==2 ? g->length_ciphertext : g->length_shared_secret; }
int shim_kem_keypair(void *k, uint8_t *pk, uint8_t *sk) { return OQS_KEM_keypair(k, pk, sk); }
int shim_kem_encaps(void *k, uint8_t *ct, uint8_t *ss, const uint8_t *pk) { return OQS_KEM_encaps(k, ct, ss, pk); }
int shim_kem_decaps(void *k, uint8_t *ss, const uint8_t *ct, const uint8_t *sk) { return OQS_KEM_decaps(k, ss, ct, sk); }
