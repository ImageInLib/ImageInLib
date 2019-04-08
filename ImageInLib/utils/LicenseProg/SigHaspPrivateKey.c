#include "stdlib.h"
#include "SigHaspPrivateKey.h"
#include <openssl/asn1.h>
#include <openssl/asn1t.h>

SIG_HASP_PKEY *SIG_HASP_PKEY_new(void)
{
	SIG_HASP_PKEY *sig;

	sig = (SIG_HASP_PKEY*)malloc(sizeof(SIG_HASP_PKEY));
	if (sig == NULL)
		return NULL;

	sig->version = 0;

	sig->p=NULL;
	sig->q=NULL;
	sig->g=NULL;

	sig->y=NULL;
	sig->x=NULL;

	return(sig);
}

void SIG_HASP_PKEY_free(SIG_HASP_PKEY *sig)
{
	if (!sig)
		return;

	if (sig->p)
		BN_clear_free(sig->p);
	if (sig->q)
		BN_clear_free(sig->q);
	if (sig->g)
		BN_clear_free(sig->g);
	if (sig->y)
		BN_clear_free(sig->y);
	if (sig->x)
		BN_clear_free(sig->x);

	free(sig);
}

static int SigHaspPrivateKey_cb(int operation, ASN1_VALUE **pval, const ASN1_ITEM *it)
{
	if(operation == ASN1_OP_NEW_PRE) {
		*pval = (ASN1_VALUE *)SIG_HASP_PKEY_new();
		if(*pval) return 2;
		return 0;
	} else if(operation == ASN1_OP_FREE_PRE) {
		SIG_HASP_PKEY_free((SIG_HASP_PKEY *)*pval);
		*pval = NULL;
		return 2;
	}
	return 1;
}

ASN1_SEQUENCE_cb(SigHaspPrivateKey, SigHaspPrivateKey_cb) = {
	ASN1_SIMPLE(SIG_HASP_PKEY, version, LONG),
	ASN1_SIMPLE(SIG_HASP_PKEY, p, BIGNUM),
	ASN1_SIMPLE(SIG_HASP_PKEY, q, BIGNUM),
	ASN1_SIMPLE(SIG_HASP_PKEY, g, BIGNUM),
	ASN1_SIMPLE(SIG_HASP_PKEY, y, BIGNUM),
	ASN1_SIMPLE(SIG_HASP_PKEY, x, BIGNUM)
} ASN1_SEQUENCE_END_cb(SIG_HASP_PKEY, SigHaspPrivateKey)

IMPLEMENT_ASN1_ENCODE_FUNCTIONS_const_fname(SIG_HASP_PKEY, SigHaspPrivateKey, SigHaspPrivateKey)
