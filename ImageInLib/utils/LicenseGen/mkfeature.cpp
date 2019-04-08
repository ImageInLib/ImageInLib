#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#ifdef _WIN32
#include "windows.h"
#else
#include "unistd.h"
typedef	unsigned long DWORD;
typedef	unsigned char BYTE;
#endif


#include <openssl/rsa.h>
#include <openssl/evp.h>
#include <openssl/objects.h>
#include <openssl/x509.h>
#include <openssl/err.h>
#include <openssl/pem.h>
#include <openssl/ssl.h>
#include <openssl/rand.h>


#define RET_OK					0
#define RET_INVALID_CMDLINE		1
#define RET_GEN_ERROR			2


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


int main(int argc, char **argv)
{
	int i, j;
	unsigned char featureID[16];
	unsigned long secretID[64];
	unsigned char secretData[64*64];
	FILE *f;

	ERR_load_crypto_strings();
	OpenSSL_add_all_algorithms();

#ifdef _WIN32
	fprintf(stdout,"Loading 'screen' into random state -");
	RAND_screen();
	fprintf(stdout," done\n");
#endif

	// generate new feature id
	if (!RAND_bytes((unsigned char*)&featureID, sizeof(featureID))) {
		fprintf(stderr,"error: Could not generate feature data");
		return RET_GEN_ERROR;
	}

	// generate new feature data
	if (!RAND_bytes((unsigned char*)&secretID, sizeof(secretID))) {
		fprintf(stderr,"error: Could not generate feature data");
		return RET_GEN_ERROR;
	}
	if (!RAND_bytes((unsigned char*)&secretData, sizeof(secretData))) {
		fprintf(stderr,"error: Could not generate feature data");
		return RET_GEN_ERROR;
	}

	f = fopen("feature.dat", "w");
	if (f == NULL) {
		fprintf(stderr,"error: Could not create feature.dat");
		return RET_GEN_ERROR;
	}

	fprintf(f, 
		"struct FEATURE {\n"
		"	char *featureName;\n"
		"	char *featureDescription;\n"
		"	unsigned char featureID[16];\n"
		"	struct {\n"
		"		unsigned long secretID;\n"
		"		unsigned char secretData[64];\n"
		"	} featureSecret[64];\n"
		"};\n"
		"\n"
	);

	fprintf(f, 
		"FEATURE _feature_ = {\n"
		"	\"FEATURE NAME\",\n"
		"	\"FEATURE DESCRIPTION\",\n"
	);
	fprintf(f, "	{");
	for (i = 0; i < 16; i++) {
		if (i != 0)
			fprintf(f, ",");
		fprintf(f, "0x%02X", featureID[i]);
	}
	fprintf(f, "},\n");
	fprintf(f, "	{\n");
	for (i = 0; i < 64; i++) {
		fprintf(f, "		{\n");
		fprintf(f, "			0x%08X,\n", secretID[i]);
		fprintf(f, "			{\n");
		for (j = 0; j < 64; j++) {
			if (j % 8 == 0)
				fprintf(f, "				");
			fprintf(f, "0x%02X", secretData[i*64+j]);
			if (j != 63)
				fprintf(f, ",");
			if (j % 8 == 7)
				fprintf(f, "\n");
		}
		fprintf(f, "			}\n");
		if (i != 63)
			fprintf(f, "		},\n");
		else
			fprintf(f, "		}\n");
	}
	fprintf(f, "	}\n");
	fprintf(f, "};\n");

	fclose(f);

	return RET_OK;	
}
