#ifndef __SIG_HASP_PKEY_H_INCLUDED__
#define __SIG_HASP_PKEY_H_INCLUDED__

#include <openssl/bn.h>

#ifdef  __cplusplus
extern "C" {
#endif

typedef struct SIG_HASP_PKEY
{
	long version;
	BIGNUM *p;
	BIGNUM *q;
	BIGNUM *g;
	BIGNUM *y;
	BIGNUM *x;
} SIG_HASP_PKEY;

SIG_HASP_PKEY * d2i_SigHaspPrivateKey(SIG_HASP_PKEY **a, const unsigned char **in, long len);
int i2d_SigHaspPrivateKey(const SIG_HASP_PKEY *a, unsigned char **out);
SIG_HASP_PKEY *SIG_HASP_PKEY_new(void);
void SIG_HASP_PKEY_free(SIG_HASP_PKEY *sig);

#ifdef  __cplusplus
}
#endif
#endif /*__SIG_HASP_PKEY_H_INCLUDED__*/
