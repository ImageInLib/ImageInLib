#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#ifdef _WIN32
#include "windows.h"
#include "io.h"
#else
#include "unistd.h"
typedef	unsigned long DWORD;
typedef	unsigned char BYTE;
#endif

#include "hasp_hl.h"

#include "SigHaspPrivateKey.h"

#include <openssl/err.h>
#include <openssl/ssl.h>
#include <openssl/rand.h>
#include <openssl/sha.h>

#include "apihook.h"
#include "crypto\crypto.h"

#include "tomocon-sig-hasp.inl"

#define RET_OK					0
#define RET_HASP_ERROR			1
#define RET_WRITE_ERROR			2
#define RET_GEN_ERROR			3
#define RET_INVALID_PEM_FILE	4


// undefine xGetSystemMetrics from MULTIMON.H
#ifdef GetSystemMetrics
#undef GetSystemMetrics
#endif

static int (WINAPI * OriginalGetSystemMetrics)(int nIndex) = GetSystemMetrics;

int WINAPI RedirectedGetSystemMetrics(int nIndex)
{
	int nRet = OriginalGetSystemMetrics(nIndex);
	if (SM_REMOTESESSION == nIndex)
		nRet = 0;
	return nRet;
};

extern "C" void hasp(int service, int seed, int lptnum, int pass1, int pass2, int *p1, int *p2, int *p3, int *p4);

bool WriteCpp(FILE *fp, const char *arrayName, BYTE *keyData, DWORD len)
{
	fprintf(fp, "unsigned char %s[] = {\n", arrayName);

	for (DWORD i = 0; i < len; i++) {
		if ((i % 12) == 0)
			fprintf(fp, "    ");
		fprintf(fp, "0x%02X", (int)keyData[i]);
		if (i == len-1) {
			fprintf(fp, "\n");
		} else {
			fprintf(fp, ",");
			if ((i % 12) == 11)
				fprintf(fp, "\n");
		}
	}

	fprintf(fp, "};\n");
	return true;
}



SIG_HASP_PKEY * ReadPrivateSigHaspKey(const char *keyfile)
{
	SIG_HASP_PKEY *pkey;
	FILE *fp;

	fp = fopen(keyfile, "r");
	if (!fp)
		return NULL;

	pkey = (SIG_HASP_PKEY *)PEM_ASN1_read((char *(*)())d2i_SigHaspPrivateKey, "HASP PRIVATE KEY", fp, NULL, 0, NULL);
	
	fclose (fp);

  	if (pkey == NULL) 
		ERR_print_errors_fp (stderr);   

	return pkey;
}


SIG_HASP_PKEY *g_tomocon_sig_hasp_PrivateKey;
typedef BYTE TCW_LICENSE_DATA[48];

bool GenerateTomoConWorkstationLicence(int nVersionNumber, DWORD dwHaspID, TCW_LICENSE_DATA *pLicenseData, SIG_HASP_PKEY *pPrivateKey)
{
	BN_MONT_CTX mont_p;
	BIGNUM r;
	BIGNUM s;
	BIGNUM kinv;
	BIGNUM k;
	BIGNUM xr;
	BIGNUM m;
	BIGNUM u1;
	BIGNUM u2;
	BIGNUM t1;
	BN_CTX *ctx = 0;
	SHA_CTX sha1;
	unsigned char buf[256];
	int iterations = 0;
	bool ret = false;

	// initialize
	BN_MONT_CTX_init(&mont_p);
	BN_init(&u1);
	BN_init(&u2);
	BN_init(&t1);
	BN_init(&kinv);
	BN_init(&k);
	BN_init(&xr);
	BN_init(&m);
	BN_init(&s);
	BN_init(&r);
	ctx = BN_CTX_new();
	if (!ctx)
		goto err;
	if (!BN_MONT_CTX_set(&mont_p, pPrivateKey->p, ctx))
		goto err;
	if (!pPrivateKey->p || !pPrivateKey->q || !pPrivateKey->g)
		goto err;
	if (BN_num_bytes(pPrivateKey->q) < 20)
		goto err;

	// sign license data
	for (;;)
	{
		for (;;)
		{
			iterations++;

			// get random k
			do {
				if (!BN_rand_range(&k, pPrivateKey->q))
					goto err;
			} while (BN_is_zero(&k));

			// compute r = (g^k mod p) mod q
			if (!BN_mod_exp_mont(&r, pPrivateKey->g, &k, pPrivateKey->p, ctx, &mont_p))
				goto err;
			if (!BN_mod(&r, &r, pPrivateKey->q, ctx))
				goto err;

			// r must be 20 bytes long
			if (BN_num_bytes(&r) != 20)
				continue;
			BN_bn2bin(&r, (*pLicenseData)+8);

			// calculate digest
			if (!CryptoGetRandom(*pLicenseData, 8))
				goto err;
			BN_bn2bin(&r, (*pLicenseData)+8);
			if (!SHA1_Init(&sha1))
				goto err;
			BN_bn2bin(pPrivateKey->y, buf);
			if (!SHA1_Update(&sha1, buf, 128))
				goto err;
			BN_bn2bin(pPrivateKey->g, buf);
			if (!SHA1_Update(&sha1, buf, 128))
				goto err;
			BN_bn2bin(pPrivateKey->p, buf);
			if (!SHA1_Update(&sha1, buf, 128))
				goto err;
			BN_bn2bin(pPrivateKey->q, buf);
			if (!SHA1_Update(&sha1, buf, 20))
				goto err;
			if (!SHA1_Update(&sha1, *pLicenseData, 28))
				goto err;
			if (!SHA1_Update(&sha1, &dwHaspID, 4))
				goto err;
			if (!SHA1_Final((*pLicenseData)+28, &sha1))
				goto err;
			memset(&sha1, 0, sizeof(sha1));
			if (BN_bin2bn((*pLicenseData)+28, 20, &m) == NULL)
				goto err;

			//  version number must match digest
			if ((*pLicenseData)[28] != nVersionNumber)
				continue;

			break;
		}

		// compute  inv(k)
		if (!BN_mod_inverse(&kinv, &k, pPrivateKey->q, ctx))
			goto err;

		// compute  s = inv(k) (m + x * r) mod q
		if (!BN_mod_mul(&xr, pPrivateKey->x,&r, pPrivateKey->q, ctx))	// s = (x * r) mod q
			goto err;	
		if (!BN_add(&s, &xr, &m))										// s = m + ((x * r) mod q)
			goto err;
		if (BN_cmp(&s, pPrivateKey->q) > 0)								// s = (m + x * r) mod q
			BN_sub(&s, &s, pPrivateKey->q);
		if (!BN_mod_mul(&s, &s, &kinv, pPrivateKey->q, ctx))			// s = inv(k) (m + x * r) mod q
			goto err;

		// s must be 20 bytes long
		if (BN_num_bytes(&s) != 20)
			continue;
		BN_bn2bin(&s, (*pLicenseData)+28);

		break;
	}

	// validate generated license data
	if (BN_bin2bn((*pLicenseData)+8, 20, &r) == NULL)
		goto err;
	if (BN_bin2bn((*pLicenseData)+28, 20, &s) == NULL)
		goto err;
	if (BN_is_zero(&r) || r.neg || BN_ucmp(&r, pPrivateKey->q) >= 0)
		goto err;
	if (BN_is_zero(&s) || s.neg || BN_ucmp(&s, pPrivateKey->q) >= 0)
		goto err;

	// t1 = inv(s)
	if ((BN_mod_inverse(&t1, &s, pPrivateKey->q, ctx)) == NULL)
		goto err;

	// u1 = (inv(s) * m) mod q
	if (!BN_mod_mul(&u1, &m, &t1, pPrivateKey->q, ctx))
		goto err;

	// u2 = (inv(s) * r) mod q
	if (!BN_mod_mul(&u2, &r, &t1, pPrivateKey->q, ctx))
		goto err;

	// t1 = ((g^u1 * y^u2) mod p) mod q
	if (!BN_mod_exp2_mont(&t1, pPrivateKey->g, &u1, pPrivateKey->y, &u2, pPrivateKey->p, ctx, &mont_p))
		goto err;
	if (!BN_mod(&t1, &t1, pPrivateKey->q, ctx))
		goto err;

	// (t1 == r) if the generated license is valid
	ret = (BN_ucmp(&t1, &r) == 0);

err:
	BN_MONT_CTX_free(&mont_p);
	BN_clear_free(&m);
	BN_clear_free(&xr);
	BN_clear_free(&k);
	BN_clear_free(&kinv);
	BN_clear_free(&r);
	BN_clear_free(&s);
	BN_clear_free(&u1);
	BN_clear_free(&u2);
	BN_clear_free(&t1);
	if (ctx)
		BN_CTX_free(ctx);
	return ret;
}

bool ValidateTomoConWorkstationLicence(DWORD dwHaspID, TCW_LICENSE_DATA *pLicenseData, int *pVersionNumber)
{
	unsigned char buf[20];
	BN_MONT_CTX mont_p;
	BIGNUM p;
	BIGNUM q;
	BIGNUM g;
	BIGNUM y;
	BIGNUM r;
	BIGNUM s;
	BIGNUM m;
	BIGNUM u1;
	BIGNUM u2;
	BIGNUM t1;
	BN_CTX *ctx = 0;
	SHA_CTX sha1;
	bool ret = false;

	// initialize
	BN_MONT_CTX_init(&mont_p);
	BN_init(&p);
	BN_init(&q);
	BN_init(&g);
	BN_init(&y);
	BN_init(&m);
	BN_init(&r);
	BN_init(&s);
	BN_init(&u1);
	BN_init(&u2);
	BN_init(&t1);
	if (BN_bin2bn(tomocon_sig_hasp_PubKey, 128, &y) == NULL)
		goto err;
	if (BN_bin2bn(tomocon_sig_hasp_PubKey + 128, 128, &g) == NULL)
		goto err;
	if (BN_bin2bn(tomocon_sig_hasp_PubKey + 2 * 128, 128, &p) == NULL)
		goto err;
	if (BN_bin2bn(tomocon_sig_hasp_PubKey + 3 * 128, 20, &q) == NULL)
		goto err;
	ctx = BN_CTX_new();
	if (!ctx)
		goto err;
	if (!BN_MONT_CTX_set(&mont_p, &p, ctx))
		goto err;

	if (BN_bin2bn((*pLicenseData)+8, 20, &r) == NULL)
		goto err;
	if (BN_is_zero(&r) || r.neg || BN_ucmp(&r, &q) >= 0)
		goto err;

	if (BN_bin2bn((*pLicenseData)+28, 20, &s) == NULL)
		goto err;
	if (BN_is_zero(&s) || s.neg || BN_ucmp(&s, &q) >= 0)
		goto err;

	SHA1_Init(&sha1);
	SHA1_Update(&sha1, tomocon_sig_hasp_PubKey, 128);
	SHA1_Update(&sha1, tomocon_sig_hasp_PubKey + 128, 128);
	SHA1_Update(&sha1, tomocon_sig_hasp_PubKey + 2 * 128, 128);
	SHA1_Update(&sha1, tomocon_sig_hasp_PubKey + 3 * 128, 20);
	SHA1_Update(&sha1, *pLicenseData, 28);
	SHA1_Update(&sha1, &dwHaspID, 4);
	SHA1_Final(buf, &sha1);
	memset(&sha1, 0, sizeof(sha1));
	if (BN_bin2bn(buf, 20, &m) == NULL)
		goto err;

	// t1 = inv(s)
	if ((BN_mod_inverse(&t1, &s, &q, ctx)) == NULL)
		goto err;

	// u1 = (inv(s) * m) mod q
	if (!BN_mod_mul(&u1, &m, &t1, &q, ctx))
		goto err;

	// u2 = (inv(s) * r) mod q
	if (!BN_mod_mul(&u2, &r, &t1, &q, ctx))
		goto err;

	// t1 = ((g^u1 * y^u2) mod p) mod q
	if (!BN_mod_exp2_mont(&t1, &g, &u1, &y, &u2, &p, ctx, &mont_p))
		goto err;
	if (!BN_mod(&t1, &t1, &q, ctx))
		goto err;

	// (t1 == r) if the generated license is valid
	ret = (BN_ucmp(&t1, &r) == 0);

	if (ret)
		*pVersionNumber = buf[0];

err:
	memset(&buf, 0, sizeof(buf));
	BN_MONT_CTX_free(&mont_p);
	BN_clear_free(&m);
	BN_clear_free(&r);
	BN_clear_free(&s);
	BN_clear_free(&u1);
	BN_clear_free(&u2);
	BN_clear_free(&t1);
	BN_clear_free(&p);
	BN_clear_free(&q);
	BN_clear_free(&g);
	BN_clear_free(&y);
	if (ctx)
		BN_CTX_free(ctx);
	return ret;
}


void NewInvalidParameterHandler(const wchar_t* expression,
   const wchar_t* function, 
   const wchar_t* file, 
   unsigned int line, 
   uintptr_t pReserved)
{
	return;
}


int main(int argc, char **argv)
{
	char buf[256];
	int p1, p2, p3, p4;
	int i;
	int nVersionNumber = 0;
	bool bHASPHL;
	bool bTimeHASP;
	int nHaspPort;
	int HaspPW1;
	int HaspPW2;
	int nHaspID;
	hasp_handle_t hHASPHL = 0;
	unsigned short HaspDATA[24];
	unsigned char *vendorCode;

	_invalid_parameter_handler oldHandler, newHandler;
	newHandler = NewInvalidParameterHandler;
	oldHandler = _set_invalid_parameter_handler(newHandler);

	if ((GetVersion() & 0x80000000) == 0) {
		if (GetSystemMetrics(SM_REMOTESESSION)) {
			CreateAPIHook(&(PVOID&)OriginalGetSystemMetrics, RedirectedGetSystemMetrics);
		}
	}

	ERR_load_crypto_strings();
	OpenSSL_add_all_algorithms();

#ifdef _WIN32
	fprintf(stdout,"Loading 'screen' into random state -");
	RAND_screen();
	fprintf(stdout," done\n");
#endif

	if (argc == 3 && _stricmp(argv[1], "-tcw") == 0) {
		nVersionNumber = atoi(argv[2]);
	}

	if (argc == 3 && _stricmp(argv[1], "-hl") == 0) {
		fprintf(stdout,"\n");
		FILE *f = fopen(argv[2], "rb");
		if (f != NULL) {
			int size = _filelength(_fileno(f));
			vendorCode = new unsigned char[size];
			if (fread(vendorCode, 1, size, f) != size) {
				delete vendorCode;
				vendorCode = NULL;
			}
			fclose(f);
		} else {
			vendorCode = NULL;
		}
		if (vendorCode == NULL) {
			fprintf(stderr,"error: Could not load HASP HL vendor code");
			return RET_HASP_ERROR;
		}
		bHASPHL = true;
	} else {
		fprintf(stdout,"\n");
		fprintf(stdout,"HASP PW1: ");
		fgets(buf, sizeof(buf), stdin);
		HaspPW1 = atoi(buf);
		fprintf(stdout,"HASP PW2: ");
		fgets(buf, sizeof(buf), stdin);
		HaspPW2 = atoi(buf);
		fprintf(stdout,"\n");
		bHASPHL = false;
	}

	if (bHASPHL) {
		// find HASP
		if (HASP_STATUS_OK != hasp_login(HASP_PROGNUM_DEFAULT_FID, vendorCode, &hHASPHL)) {
			fprintf(stderr, "error: HASP not found");
			return RET_HASP_ERROR;
		}
		char *info = 0;
		if (HASP_STATUS_OK !=hasp_get_sessioninfo(hHASPHL, HASP_KEYINFO, &info)) {
			fprintf(stderr, "error: HASP not found");
			return RET_HASP_ERROR;
		}
		if (strstr(info, "<rtc/>") != NULL)
			bTimeHASP = true;
		else
			bTimeHASP = false;
		// get HASPID
		char *p;
		p = strstr(info, "<haspid>");
		if (p == NULL) {
			hasp_free(info);
			fprintf(stderr, "error: HASP not found");
			return RET_HASP_ERROR;
		}
		p += strlen("<haspid>");
		nHaspID = atol(p);
		hasp_free(info);
	} else {
		// find HASP
		p1 = p2 = p3 = p4 = 0;
		hasp(5, 0, 0, HaspPW1, HaspPW2, &p1, &p2, &p3, &p4);
		if ((p1 != 4 || p2 != 1) && (p1 != 0 || p2 != 5)) {
			fprintf(stderr, "error: HASP not found");
			return RET_HASP_ERROR;
		}
		nHaspPort = p3;
		if (p2 == 5)
			bTimeHASP = true;
		else
			bTimeHASP = false;
		// get HASPID
		p1 = p2 = p3 = p4 = 0;
		hasp(6, 0, nHaspPort, HaspPW1, HaspPW2, &p1, &p2, &p3, &p4);
		if (p3 != 0) {
			fprintf(stderr, "error: HASP not found");
			return RET_HASP_ERROR;
		}
		nHaspID = p1 + (p2 << 16);
	}

	// get user confirmation
	if (bHASPHL)
		fprintf(stdout,"HASP HL%s [ID=%08X] found.\n", bTimeHASP ? " Time" : "", nHaspID);
	else
		fprintf(stdout,"%sHASP M4 [ID=%08X] found.\n", bTimeHASP ? "Time" : "", nHaspID);
	fprintf(stdout,"Are you sure, you want to overwrite HASP memory? [Y/N] ");
	fgets(buf, sizeof(buf), stdin);
	if (strcmp(buf, "y\n") != 0 && strcmp(buf, "Y\n") != 0)
		return RET_OK;
	fprintf(stdout,"\n");

	// generate new memory contents
	if (nVersionNumber > 0 && nVersionNumber <= 32)
	{
		fprintf(stdout,"Loading tomocon-sig-hasp private key\n");
		g_tomocon_sig_hasp_PrivateKey = ReadPrivateSigHaspKey("tomocon-sig-hasp.pri");
		if (g_tomocon_sig_hasp_PrivateKey == NULL)
			return RET_INVALID_PEM_FILE;

		if (!GenerateTomoConWorkstationLicence(nVersionNumber, nHaspID, (TCW_LICENSE_DATA*)&HaspDATA, g_tomocon_sig_hasp_PrivateKey))
		{
			fprintf(stderr,"error: Could not generate TomoCon Workstation license data");
			return RET_GEN_ERROR;
		}

		int nVersion;
		if (!ValidateTomoConWorkstationLicence(nHaspID, (TCW_LICENSE_DATA*)&HaspDATA, &nVersion) || nVersion != nVersionNumber)
		{
			fprintf(stderr,"error: Could not generate TomoCon Workstation license data");
			return RET_GEN_ERROR;
		}
	}
	else
	{
		if (!RAND_bytes((unsigned char*)&HaspDATA, sizeof(HaspDATA)))
		{
			fprintf(stderr,"error: Could not generate license data");
			return RET_GEN_ERROR;
		}
	}

	// write to HASP
	fprintf(stdout,"Writing to HASP memory -");
	if (bHASPHL) {
		if (HASP_STATUS_OK != hasp_write(hHASPHL, HASP_FILEID_MAIN, 0, 48, HaspDATA)) {
			fprintf(stdout," FAILED\n");
			fprintf(stderr,"error: could not write to HASP");
			return RET_HASP_ERROR;
		}
	} else {
		p1 = 0;
		p2 = 24;
		p3 = 0;
		p4 = (int)&HaspDATA;
		hasp(51, 0, nHaspPort, HaspPW1, HaspPW2, &p1, &p2, &p3, &p4);
		if (p3 != 0) {
			fprintf(stdout," FAILED\n");
			fprintf(stderr,"error: could not write to HASP");
			return RET_HASP_ERROR;
		}
	}
	fprintf(stdout," OK\n");

	// write data file
	FILE *f;
	sprintf(buf, "%08X.HSP", nHaspID);
	f = fopen(buf, "w");
	if (f == NULL) {
		fprintf(stderr,"error: could not create file %s", buf);
		return RET_WRITE_ERROR;
	}
	fprintf(f, "HASPID = \"%08X\"\n", nHaspID);
	fprintf(f, "HASPData = \"");
	for (i = 0; i < 24; i++) {
		if (i != 0)
			fprintf(f, " ");
		fprintf(f, "0x%04x", (int)HaspDATA[i]);
	}
	fprintf(f, "\"\n");
	if (nVersionNumber > 0 && nVersionNumber <= 32)
		fprintf(f, "VersionNumber = \"%d\"", nVersionNumber);
	fclose(f);

	return RET_OK;	
}
