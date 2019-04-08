#include "headers.h"
#include "ini.h"
#include "TaggedDataBuffer.h"
#include "License.h"
#include "builtinkeys.h"
#include "datetime.h"
#include "mutex.h"
#include "PRF.h"

#include "hasp_hl.h"

#include <openssl/rsa.h>
#include <openssl/evp.h>
#include <openssl/objects.h>
#include <openssl/x509.h>
#include <openssl/err.h>
#include <openssl/pem.h>
#include <openssl/ssl.h>
#include <openssl/bio.h>
#include <openssl/rand.h>
#include "apihook.h"

#include "features.h"

#define RET_OK					0
#define RET_INVALID_CMDLINE		1
#define RET_INVALID_PEM_FILE	3
#define RET_WRITE_ERROR			4
#define RET_INVALID_DEF_FILE	5
#define RET_MISSING_HASP		6

static int (WINAPI * OriginalGetSystemMetrics)(int nIndex) = GetSystemMetrics;

int WINAPI RedirectedGetSystemMetrics(int nIndex)
{
	int nRet = OriginalGetSystemMetrics(nIndex);
	if (SM_REMOTESESSION == nIndex)
		nRet = 0;
	return nRet;
};

extern "C" void hasp(int service, int seed, int lptnum, int pass1, int pass2, int *p1, int *p2, int *p3, int *p4);

const char *g_inputFile;
char g_outputFile[1024];


RSA *g_tomocon_lic_rsa;
DSA *g_tomocon_sig_lic_dsa;
DSA *g_tomocon_sig_ini_dsa;

RSA *g_pacs_lic_rsa;
DSA *g_pacs_sig_lic_dsa;
DSA *g_pacs_sig_ini_dsa;

enum LicenceApplicationType {
	TomoConWorkstation,
	PACSServer,
	miniPACSServer
};
LicenceApplicationType g_application_type;


const char *g_IssuedBy = "TatraMed s.r.o.";
const char *g_IssuedTo = NULL;
const char *g_ValidFrom = NULL;
const char *g_ValidTo = NULL;
const char *g_HASPID = NULL;
const char *g_Comment = NULL;
const FEATURE *g_Feature = NULL;
BYTE *g_FeatureData = NULL;
DWORD g_FeatureDataSize;


int g_HASPPW1 = 0;
int g_HASPPW2 = 0;

unsigned char *g_HASPHL_vendor_code = NULL;

DWORD g_dwHASPID;

unsigned short g_HASPData[24];
unsigned short g_HASPCrypt[24];


static char *usage_text = 
	"Usage:\n"
	"  LicenseGen license.def\n"
	"\n"
	"";

static void usage() 
{
	printf(usage_text);
}

DSA * ReadDSAKey(const char *keyfile, bool pub)
{
	DSA *pkey;
	FILE *fp;

	fp = fopen(keyfile, "r");
	if (!fp)
		return NULL;

	if (pub) {

		EVP_PKEY *pk;
		X509 *x509 = PEM_read_X509(fp, NULL, 0, NULL);		
		fclose (fp);
		if (x509 == NULL) {  
			ERR_print_errors_fp (stderr);
			return NULL;   
		}
		pk = X509_get_pubkey(x509);
		X509_free(x509);

		pkey = EVP_PKEY_get1_DSA(pk);
		
	} else {

		pkey = (DSA *)PEM_ASN1_read((char *(*)())(pub ? d2i_DSAPublicKey : d2i_DSAPrivateKey), pub ? PEM_STRING_PUBLIC : PEM_STRING_DSA, fp, NULL, 0, NULL);
		fclose (fp);

	}

  	if (pkey == NULL) 
		ERR_print_errors_fp (stderr);   

	return pkey;
}


RSA * ReadRSAKey(const char *keyfile, bool pub)
{
	RSA *pkey;
	FILE *fp;

	fp = fopen(keyfile, "r");
	if (!fp)
		return NULL;

	if (pub) {

		EVP_PKEY *pk;
		X509 *x509 = PEM_read_X509(fp, NULL, 0, NULL);		
		fclose (fp);
		if (x509 == NULL) {  
			ERR_print_errors_fp (stderr);
			return NULL;   
		}
		pk = X509_get_pubkey(x509);
		X509_free(x509);

		pkey = EVP_PKEY_get1_RSA(pk);
		
	} else {

		pkey = (RSA *)PEM_ASN1_read((char *(*)())(pub ? d2i_RSAPublicKey : d2i_RSAPrivateKey), pub ? PEM_STRING_PUBLIC : PEM_STRING_RSA, fp, NULL, 0, NULL);
		fclose (fp);

	}

  	if (pkey == NULL) 
		ERR_print_errors_fp (stderr);   

	return pkey;
}


INISTATUS SignLicense(
	VOID       *pCallbackCtx,
	DWORD        dwDigestType,
	const BYTE  *pDigest,
	DWORD		 cbDigest,
	DWORD       *pdwSigType,
	BYTE       **ppPubKeyID,	// must be allocated by malloc !
	DWORD		*pcbPubKeyID,
	BYTE       **ppSignature,   // must be allocated by malloc !
	DWORD		*pcbSignature
	)
{
	static BYTE buf[16384];
	DWORD len;

	if (g_application_type == TomoConWorkstation) {
		*pdwSigType = INI_SIGNATURE_DSA;
		*ppPubKeyID = (BYTE*)malloc(tomocon_sig_ini_dsa_KeyID_size);
		*pcbPubKeyID = tomocon_sig_ini_dsa_KeyID_size;
		memmove(*ppPubKeyID, tomocon_sig_ini_dsa_KeyID, tomocon_sig_ini_dsa_KeyID_size);
	} else if (g_application_type == PACSServer || g_application_type == miniPACSServer) {
		*pdwSigType = INI_SIGNATURE_DSA;
		*ppPubKeyID = (BYTE*)malloc(tm_pacs_sig_ini_dsa_KeyID_size);
		*pcbPubKeyID = tm_pacs_sig_ini_dsa_KeyID_size;
		memmove(*ppPubKeyID, tm_pacs_sig_ini_dsa_KeyID, tm_pacs_sig_ini_dsa_KeyID_size);
	}

	HCRYPTOSIGN hSig;
	DWORD dw;
	CryptoCreateSignature(CRYPTO_ALG_DSA, &hSig);

	if (g_application_type == TomoConWorkstation) {
		len = BN_num_bytes(g_tomocon_sig_ini_dsa->p);
		BN_bn2bin(g_tomocon_sig_ini_dsa->p, buf);
		CryptoSetSignatureParam(hSig, CRYPTO_SP_DSA_P, buf, len);	
		len = BN_num_bytes(g_tomocon_sig_ini_dsa->q);
		BN_bn2bin(g_tomocon_sig_ini_dsa->q, buf);
		CryptoSetSignatureParam(hSig, CRYPTO_SP_DSA_Q, buf, len);
		len = BN_num_bytes(g_tomocon_sig_ini_dsa->g);
		BN_bn2bin(g_tomocon_sig_ini_dsa->g, buf);
		CryptoSetSignatureParam(hSig, CRYPTO_SP_DSA_G, buf, len);
		len = BN_num_bytes(g_tomocon_sig_ini_dsa->priv_key);
		BN_bn2bin(g_tomocon_sig_ini_dsa->priv_key, buf);
		CryptoSetSignatureParam(hSig, CRYPTO_SP_DSA_PRIV, buf, len);
		len = BN_num_bytes(g_tomocon_sig_ini_dsa->pub_key);
		BN_bn2bin(g_tomocon_sig_ini_dsa->pub_key, buf);
		CryptoSetSignatureParam(hSig, CRYPTO_SP_DSA_PUB, buf, len);
	} else if (g_application_type == PACSServer || g_application_type == miniPACSServer) {
		len = BN_num_bytes(g_pacs_sig_ini_dsa->p);
		BN_bn2bin(g_pacs_sig_ini_dsa->p, buf);
		CryptoSetSignatureParam(hSig, CRYPTO_SP_DSA_P, buf, len);	
		len = BN_num_bytes(g_pacs_sig_ini_dsa->q);
		BN_bn2bin(g_pacs_sig_ini_dsa->q, buf);
		CryptoSetSignatureParam(hSig, CRYPTO_SP_DSA_Q, buf, len);
		len = BN_num_bytes(g_pacs_sig_ini_dsa->g);
		BN_bn2bin(g_pacs_sig_ini_dsa->g, buf);
		CryptoSetSignatureParam(hSig, CRYPTO_SP_DSA_G, buf, len);
		len = BN_num_bytes(g_pacs_sig_ini_dsa->priv_key);
		BN_bn2bin(g_pacs_sig_ini_dsa->priv_key, buf);
		CryptoSetSignatureParam(hSig, CRYPTO_SP_DSA_PRIV, buf, len);
		len = BN_num_bytes(g_pacs_sig_ini_dsa->pub_key);
		BN_bn2bin(g_pacs_sig_ini_dsa->pub_key, buf);
		CryptoSetSignatureParam(hSig, CRYPTO_SP_DSA_PUB, buf, len);
	}

	if (dwDigestType == INI_CHECKSUM_MD5)
		dw = CRYPTO_DIGEST_MD5;
	else if (dwDigestType == INI_CHECKSUM_SHA1)
		dw = CRYPTO_DIGEST_SHA1;
	else {
		dw = CRYPTO_DIGEST_SHA1;
		pDigest += 16;
		cbDigest = 20;
	}

	CryptoSetSignatureParam(hSig, CRYPTO_SP_DIGESTTYPE, (BYTE*)&dw, sizeof(dw));
	CryptoSetSignatureParam(hSig, CRYPTO_SP_DIGESTVALUE, pDigest, cbDigest);
	*ppSignature = (BYTE*)malloc(1024);
	*pcbSignature = 1024;
	CryptoGetSignatureParam(hSig, CRYPTO_SP_SIGNATURE, *ppSignature, pcbSignature);
	CryptoDestroySignature(hSig);

	return STATUS_INI_SUCCESS;
}

void binDump(const char *pszName, const BYTE* pData, DWORD dwLength)
{
	printf("\n%s:", pszName);
	for (DWORD i = 0; i < dwLength; i++) {
		if ((i % 32) == 0)
			printf("\n");
		printf("%02X", pData[i]);
		if ((i % 32) != 31)
			printf(" ");
	}
}


#include <crtdbg.h>

void CreateCertificate()
{
	static BYTE buf[16384];
	HTDBELEMENT hMainSeq;
	HTDBITEM hMainItem;
	HTDBELEMENT hHaspSeq;
	HTDBITEM hHaspItem;
	HTDBELEMENT hSignSeq;
	HTDBITEM hSignItem;
	CTaggedDataBuffer lic;
	BYTE licID_bin[20];
	char licID_b64[40];
	HANDLE hLicFile;
	DWORD dwLength;
	BYTE *pData;
	DWORD dwDataLength;
	char *pDataBase64;
	DWORD dwDataBase64Length;
	DWORD dwVersion;
	DWORD dwType;
	DWORD dwAlgo;
	HCRYPTOSIGN hSig;
	Security::CBinary encryptedLicense;
	Security::CBinary licenseBlock;
	Security::CBinary rsaBlock;
	Security::CBinary haspSalt;
	DWORD dwValue;
	BYTE AES_KEY[16+16];
	BYTE AES_KEY_256[32+16];
	HCRYPTOCIPHER hAES;
	HCRYPTOCIPHER hRSA;
	INISTATUS status;
	int pad;
	int i;

	CryptoGetRandom(licID_bin, sizeof(licID_bin));
	dwLength = sizeof(licID_b64);
	CryptoConvertToBase64(licID_bin, sizeof(licID_bin), licID_b64, &dwLength, false);

	// create LICENSE_TAG_MAIN_SEQ
	hMainSeq = lic.CreateSequence(NULL, LICENSE_TAG_MAIN_SEQ);
	hMainItem = lic.CreateSequenceItem(hMainSeq);

	// set LICENSE_TAG_MAIN_SEQ
	dwVersion = 1;
	lic.SetElementValue(hMainItem, CTaggedDataBuffer::etValue, LICENSE_TAG_VERSION, sizeof(dwVersion), (BYTE*)&dwVersion);
	dwType = LICENSE_TYPE_HASP;
	lic.SetElementValue(hMainItem, CTaggedDataBuffer::etValue, LICENSE_TAG_TYPE, sizeof(DWORD), (BYTE*)&dwType);
	lic.SetElementValue(hMainItem, CTaggedDataBuffer::etValue, LICENSE_TAG_ID, sizeof(licID_bin), licID_bin);
	if (g_application_type == TomoConWorkstation)
		lic.SetElementValue(hMainItem, CTaggedDataBuffer::etValue, LICENSE_TAG_APPLICATION, tomocon_AppID_size, tomocon_AppID);
	else if (g_application_type == PACSServer)
		lic.SetElementValue(hMainItem, CTaggedDataBuffer::etValue, LICENSE_TAG_APPLICATION, tm_pacs_AppID_size, tm_pacs_AppID);
	else if (g_application_type == miniPACSServer)
		lic.SetElementValue(hMainItem, CTaggedDataBuffer::etValue, LICENSE_TAG_APPLICATION, tm_minipacs_AppID_size, tm_minipacs_AppID);
	lic.SetElementValue(hMainItem, CTaggedDataBuffer::etValue, LICENSE_TAG_ISSUER, strlen(g_IssuedBy), (BYTE*)g_IssuedBy);
	lic.SetElementValue(hMainItem, CTaggedDataBuffer::etValue, LICENSE_TAG_FEATURE, strlen(g_Feature->featureName), (BYTE*)g_Feature->featureName);
	lic.SetElementValue(hMainItem, CTaggedDataBuffer::etValue, LICENSE_TAG_FEATURE_DESC, strlen(g_Feature->featureDescription), (BYTE*)g_Feature->featureDescription);
	if (g_FeatureData && g_FeatureDataSize)
		lic.SetElementValue(hMainItem, CTaggedDataBuffer::etValue, LICENSE_TAG_FEATURE_DATA, g_FeatureDataSize, g_FeatureData);
	lic.SetElementValue(hMainItem, CTaggedDataBuffer::etValue, LICENSE_TAG_FEATURE_ID, sizeof(g_Feature->featureID), (BYTE*)g_Feature->featureID);
	if (g_IssuedTo)
		lic.SetElementValue(hMainItem, CTaggedDataBuffer::etValue, LICENSE_TAG_SUBJECT, strlen(g_IssuedTo), (BYTE*)g_IssuedTo);
	if (g_Comment)
		lic.SetElementValue(hMainItem, CTaggedDataBuffer::etValue, LICENSE_TAG_COMMENT, strlen(g_Comment), (BYTE*)g_Comment);
	if (g_ValidFrom) {
		long date = ConvertStringToDate(g_ValidFrom);
		lic.SetElementValue(hMainItem, CTaggedDataBuffer::etValue, LICENSE_TAG_VALID_FROM, sizeof(date), (BYTE*)&date);
	}
	if (g_ValidTo) {
		long date = ConvertStringToDate(g_ValidTo);
		lic.SetElementValue(hMainItem, CTaggedDataBuffer::etValue, LICENSE_TAG_VALID_TO, sizeof(date), (BYTE*)&date);
	}

	// create LICENSE_TAG_HASP_SEQ
	hHaspSeq = lic.CreateSequence(hMainItem, LICENSE_TAG_HASP_SEQ);
	hHaspItem = lic.CreateSequenceItem(hHaspSeq);

	// set LICENSE_TAG_HASP_SEQ
	lic.SetElementValue(hHaspItem, CTaggedDataBuffer::etValue, LICENSE_TAG_HASP_HASPID, sizeof(g_dwHASPID), (BYTE*)&g_dwHASPID);
	haspSalt.resize(64);
	CryptoGetRandom((BYTE*)(const BYTE*)haspSalt, haspSalt.size());
	lic.SetElementValue(hHaspItem, CTaggedDataBuffer::etValue, LICENSE_TAG_HASP_SALT, haspSalt.size(), (BYTE*)(const BYTE*)haspSalt);

	// set LICENSE_TAG_LICENSE_BLOCK
	Security::CBinary secretData;
	secretData.resize(64 * (sizeof(g_Feature->featureSecret[0].secretID) + sizeof(g_Feature->featureSecret[0].secretData)));
	for (i = 0; i < 64; i++) {
		BYTE *p = ((BYTE*)(const BYTE*)secretData) + i * (4 + 64);
		memmove(p, &g_Feature->featureSecret[i].secretID, 4);
		p += 4;
		memmove(p, &g_Feature->featureSecret[i].secretData, 64);
		HCRYPTOCIPHER hAES;
		CryptoCreateCipher(CRYPTO_ALG_AES, &hAES);
		dwValue = CRYPTO_MODE_CBC;
		CryptoSetCipherParam(hAES, CRYPTO_CP_MODE, (BYTE*)&dwValue, sizeof(DWORD));
		Security::PRF_Expand((const BYTE*)g_HASPData, sizeof(g_HASPData), 
			(const BYTE*)&g_Feature->featureSecret[i].secretID, sizeof(g_Feature->featureSecret[i].secretID), 
			(const BYTE*)haspSalt, haspSalt.size(), 
			AES_KEY, sizeof(AES_KEY));
		CryptoSetCipherParam(hAES, CRYPTO_CP_KEY, (BYTE*)AES_KEY, 16);
		CryptoSetCipherParam(hAES, CRYPTO_CP_IV, (BYTE*)AES_KEY+16, 16);
		CryptoEncrypt(hAES, NULL, p, 64);
		CryptoDestroyCipher(hAES);
	}
	lic.SetElementValue(hMainItem, CTaggedDataBuffer::etValue, LICENSE_TAG_SECRET_BLOCK, secretData.size(), (const BYTE*)secretData);

	// create CERT_TAG_SIGNATURE_SEQ
	hSignSeq = lic.CreateSequence(NULL, LICENSE_TAG_SIGNATURE_SEQ);
	hSignItem = lic.CreateSequenceItem(hSignSeq);

	// sign CERT_TAG_MAIN_SEQ
	CryptoCreateSignature(CRYPTO_ALG_DSA, &hSig);
	dwAlgo = CRYPTO_DIGEST_SHA1;
	CryptoSetSignatureParam(hSig, CRYPTO_SP_DIGESTTYPE, (BYTE*)&dwAlgo, sizeof(dwAlgo));
	if (g_application_type == TomoConWorkstation) {
		dwLength = BN_num_bytes(g_tomocon_sig_lic_dsa->p);
		BN_bn2bin(g_tomocon_sig_lic_dsa->p, buf);
		CryptoSetSignatureParam(hSig, CRYPTO_SP_DSA_P, buf, dwLength);	
		dwLength = BN_num_bytes(g_tomocon_sig_lic_dsa->q);
		BN_bn2bin(g_tomocon_sig_lic_dsa->q, buf);
		CryptoSetSignatureParam(hSig, CRYPTO_SP_DSA_Q, buf, dwLength);
		dwLength = BN_num_bytes(g_tomocon_sig_lic_dsa->g);
		BN_bn2bin(g_tomocon_sig_lic_dsa->g, buf);
		CryptoSetSignatureParam(hSig, CRYPTO_SP_DSA_G, buf, dwLength);
		dwLength = BN_num_bytes(g_tomocon_sig_lic_dsa->priv_key);
		BN_bn2bin(g_tomocon_sig_lic_dsa->priv_key, buf);
		CryptoSetSignatureParam(hSig, CRYPTO_SP_DSA_PRIV, buf, dwLength);
		dwLength = BN_num_bytes(g_tomocon_sig_lic_dsa->pub_key);
		BN_bn2bin(g_tomocon_sig_lic_dsa->pub_key, buf);
		CryptoSetSignatureParam(hSig, CRYPTO_SP_DSA_PUB, buf, dwLength);
	} else if (g_application_type == PACSServer || g_application_type == miniPACSServer) {
		dwLength = BN_num_bytes(g_pacs_sig_lic_dsa->p);
		BN_bn2bin(g_pacs_sig_lic_dsa->p, buf);
		CryptoSetSignatureParam(hSig, CRYPTO_SP_DSA_P, buf, dwLength);	
		dwLength = BN_num_bytes(g_pacs_sig_lic_dsa->q);
		BN_bn2bin(g_pacs_sig_lic_dsa->q, buf);
		CryptoSetSignatureParam(hSig, CRYPTO_SP_DSA_Q, buf, dwLength);
		dwLength = BN_num_bytes(g_pacs_sig_lic_dsa->g);
		BN_bn2bin(g_pacs_sig_lic_dsa->g, buf);
		CryptoSetSignatureParam(hSig, CRYPTO_SP_DSA_G, buf, dwLength);
		dwLength = BN_num_bytes(g_pacs_sig_lic_dsa->priv_key);
		BN_bn2bin(g_pacs_sig_lic_dsa->priv_key, buf);
		CryptoSetSignatureParam(hSig, CRYPTO_SP_DSA_PRIV, buf, dwLength);
		dwLength = BN_num_bytes(g_pacs_sig_lic_dsa->pub_key);
		BN_bn2bin(g_pacs_sig_lic_dsa->pub_key, buf);
		CryptoSetSignatureParam(hSig, CRYPTO_SP_DSA_PUB, buf, dwLength);
	}
	lic.GetElementInfo(hMainSeq, NULL, NULL, &dwLength, &pData);
	CryptoSignatureHashData(hSig, pData, dwLength);
	dwLength = sizeof(buf);
	pData = buf;
	CryptoGetSignatureParam(hSig, CRYPTO_SP_SIGNATURE, pData, &dwLength);
	CryptoDestroySignature(hSig);

	// set LICENSE_TAG_SIGNATURE_SEQ item
	dwAlgo = LICENSE_SIGNATURE_DSA_SHA1;
	lic.SetElementValue(hSignItem, CTaggedDataBuffer::etValue, LICENSE_TAG_SIGN_ALGO, sizeof(dwAlgo), (BYTE*)&dwAlgo);
	if (g_application_type == TomoConWorkstation)
		lic.SetElementValue(hSignItem, CTaggedDataBuffer::etValue, LICENSE_TAG_SIGN_KEY_ID, tomocon_sig_lic_dsa_KeyID_size, tomocon_sig_lic_dsa_KeyID);
	else if (g_application_type == PACSServer || g_application_type == miniPACSServer)
		lic.SetElementValue(hSignItem, CTaggedDataBuffer::etValue, LICENSE_TAG_SIGN_KEY_ID, tm_pacs_sig_lic_dsa_KeyID_size, tm_pacs_sig_lic_dsa_KeyID);
	pData = buf;
	lic.SetElementValue(hSignItem, CTaggedDataBuffer::etValue, LICENSE_TAG_SIGN_DATA, dwLength, pData);

	// create certificate file
	IniCreateFile(&hLicFile);

	// write certificate info
	IniSetValue(hLicFile, "Version", INI_SZ, "1");
	IniSetValue(hLicFile, "LicenseID", INI_SZ, licID_b64);
	IniSetValue(hLicFile, "LicenseType", INI_SZ, "HASP");
	IniSetValue(hLicFile, "HASPID", INI_SZ, g_HASPID);
	if (g_application_type == TomoConWorkstation)
		IniSetValue(hLicFile, "Application", INI_SZ, "TomoCon 3.0 Workstation");
	else if (g_application_type == PACSServer)
		IniSetValue(hLicFile, "Application", INI_SZ, "TomoCon PACS Server");
	else if (g_application_type == miniPACSServer)
		IniSetValue(hLicFile, "Application", INI_SZ, "TomoCon miniPACS Server");
	IniSetValue(hLicFile, "IssedBy", INI_SZ, "TatraMed s.r.o.");
	if (g_IssuedTo)
		IniSetValue(hLicFile, "IssedTo", INI_SZ, g_IssuedTo);
	if (g_ValidFrom)
		IniSetValue(hLicFile, "ValidFrom", INI_SZ, g_ValidFrom);
	if (g_ValidTo)
		IniSetValue(hLicFile, "ValidTo", INI_SZ, g_ValidTo);
	IniSetValue(hLicFile, "Feature", INI_SZ, g_Feature->featureName);
	IniSetValue(hLicFile, "FeatureDescription", INI_SZ, g_Feature->featureDescription);
	if (g_Comment)
		IniSetValue(hLicFile, "Comment", INI_SZ, g_Comment);

	// write certificate data
	lic.Write(NULL, &dwDataLength);
	pData = new BYTE[dwDataLength];
	lic.Write(pData, &dwDataLength);

	// encrypt license block
	dwLength = dwDataLength;
	licenseBlock.resize(48 + 16 + (dwLength & ~15) + 16);
	CryptoGetRandom(licenseBlock, 48 + 16);
	memmove(licenseBlock + 48 + 16, pData, dwLength);
	delete pData;
	pad = (dwLength & ~15) + 16 - dwLength;
	memset(licenseBlock + 48 + 16 + dwLength, pad, pad); 
	CryptoCreateCipher(CRYPTO_ALG_AES, &hAES);
	dwValue = CRYPTO_MODE_CBC;
	CryptoSetCipherParam(hAES, CRYPTO_CP_MODE, (BYTE*)&dwValue, sizeof(DWORD));
	Security::PRF_Expand((const BYTE*)g_HASPData, sizeof(g_HASPData),
		(const BYTE*)&g_dwHASPID, sizeof(g_dwHASPID),
		licenseBlock, 48, 
		AES_KEY_256, sizeof(AES_KEY_256));
	dwValue = 32;
	CryptoSetCipherParam(hAES, CRYPTO_CP_KEYLEN, (BYTE*)&dwValue, sizeof(DWORD));
	CryptoSetCipherParam(hAES, CRYPTO_CP_KEY, (BYTE*)AES_KEY_256, 32);
	CryptoSetCipherParam(hAES, CRYPTO_CP_IV, (BYTE*)AES_KEY_256+32, 16);
	CryptoEncrypt(hAES, NULL, licenseBlock + 48, licenseBlock.size() - 48);
	CryptoDestroyCipher(hAES);

	// RSA encrypt license block
	for (i = 0; i < 48 / sizeof(unsigned short); i++) {
		register unsigned short *p = (unsigned short*)((const BYTE*)licenseBlock) + i;
		*p ^= g_HASPData[i];
	}
	for (i = 0; i < 48 / sizeof(DWORD); i++) {
		register DWORD *p = (DWORD*)((const BYTE*)licenseBlock) + i;
		*p ^= g_dwHASPID;
	}
	CryptoCreateCipher(CRYPTO_ALG_RSA, &hRSA);
	dwValue = CRYPTO_MODE_PRIV;
	CryptoSetCipherParam(hRSA, CRYPTO_CP_MODE, (BYTE*)&dwValue, sizeof(DWORD));	
	if (g_application_type == TomoConWorkstation) {
		dwLength = BN_num_bytes(g_tomocon_lic_rsa->n);
		BN_bn2bin(g_tomocon_lic_rsa->n, buf);
		CryptoSetCipherParam(hRSA, CRYPTO_CP_RSA_N, buf, dwLength);	
		dwLength = BN_num_bytes(g_tomocon_lic_rsa->e);
		BN_bn2bin(g_tomocon_lic_rsa->e, buf);
		CryptoSetCipherParam(hRSA, CRYPTO_CP_RSA_E, buf, dwLength);	
		dwLength = BN_num_bytes(g_tomocon_lic_rsa->d);
		BN_bn2bin(g_tomocon_lic_rsa->d, buf);
		CryptoSetCipherParam(hRSA, CRYPTO_CP_RSA_D, buf, dwLength);	
		dwLength = BN_num_bytes(g_tomocon_lic_rsa->p);
		BN_bn2bin(g_tomocon_lic_rsa->p, buf);
		CryptoSetCipherParam(hRSA, CRYPTO_CP_RSA_P, buf, dwLength);	
		dwLength = BN_num_bytes(g_tomocon_lic_rsa->q);
		BN_bn2bin(g_tomocon_lic_rsa->q, buf);
		CryptoSetCipherParam(hRSA, CRYPTO_CP_RSA_Q, buf, dwLength);	
		if (g_tomocon_lic_rsa->dmp1) {
			dwLength = BN_num_bytes(g_tomocon_lic_rsa->dmp1);
			BN_bn2bin(g_tomocon_lic_rsa->dmp1, buf);
			CryptoSetCipherParam(hRSA, CRYPTO_CP_RSA_DMP1, buf, dwLength);	
		}
		if (g_tomocon_lic_rsa->dmq1) {
			dwLength = BN_num_bytes(g_tomocon_lic_rsa->dmq1);
			BN_bn2bin(g_tomocon_lic_rsa->dmq1, buf);
			CryptoSetCipherParam(hRSA, CRYPTO_CP_RSA_DMQ1, buf, dwLength);	
		}
		if (g_tomocon_lic_rsa->iqmp) {
			dwLength = BN_num_bytes(g_tomocon_lic_rsa->iqmp);
			BN_bn2bin(g_tomocon_lic_rsa->iqmp, buf);
			CryptoSetCipherParam(hRSA, CRYPTO_CP_RSA_IQMP, buf, dwLength);	
		}
	} else if (g_application_type == PACSServer || g_application_type == miniPACSServer) {
		dwLength = BN_num_bytes(g_pacs_lic_rsa->n);
		BN_bn2bin(g_pacs_lic_rsa->n, buf);
		CryptoSetCipherParam(hRSA, CRYPTO_CP_RSA_N, buf, dwLength);	
		dwLength = BN_num_bytes(g_pacs_lic_rsa->e);
		BN_bn2bin(g_pacs_lic_rsa->e, buf);
		CryptoSetCipherParam(hRSA, CRYPTO_CP_RSA_E, buf, dwLength);	
		dwLength = BN_num_bytes(g_pacs_lic_rsa->d);
		BN_bn2bin(g_pacs_lic_rsa->d, buf);
		CryptoSetCipherParam(hRSA, CRYPTO_CP_RSA_D, buf, dwLength);	
		dwLength = BN_num_bytes(g_pacs_lic_rsa->p);
		BN_bn2bin(g_pacs_lic_rsa->p, buf);
		CryptoSetCipherParam(hRSA, CRYPTO_CP_RSA_P, buf, dwLength);	
		dwLength = BN_num_bytes(g_pacs_lic_rsa->q);
		BN_bn2bin(g_pacs_lic_rsa->q, buf);
		CryptoSetCipherParam(hRSA, CRYPTO_CP_RSA_Q, buf, dwLength);	
		if (g_pacs_lic_rsa->dmp1) {
			dwLength = BN_num_bytes(g_pacs_lic_rsa->dmp1);
			BN_bn2bin(g_pacs_lic_rsa->dmp1, buf);
			CryptoSetCipherParam(hRSA, CRYPTO_CP_RSA_DMP1, buf, dwLength);	
		}
		if (g_pacs_lic_rsa->dmq1) {
			dwLength = BN_num_bytes(g_pacs_lic_rsa->dmq1);
			BN_bn2bin(g_pacs_lic_rsa->dmq1, buf);
			CryptoSetCipherParam(hRSA, CRYPTO_CP_RSA_DMQ1, buf, dwLength);	
		}
		if (g_pacs_lic_rsa->iqmp) {
			dwLength = BN_num_bytes(g_pacs_lic_rsa->iqmp);
			BN_bn2bin(g_pacs_lic_rsa->iqmp, buf);
			CryptoSetCipherParam(hRSA, CRYPTO_CP_RSA_IQMP, buf, dwLength);	
		}
	}
	dwValue = CRYPTO_RSA_PKCS1_PADDING;
	CryptoSetCipherParam(hRSA, CRYPTO_CP_RSA_PADDING, (BYTE*)&dwValue, sizeof(DWORD));	
	rsaBlock.resize(16384);
	dwLength = 16384;
	CryptoEncrypt(hRSA, NULL, licenseBlock, 512-16, rsaBlock, &dwLength);
	rsaBlock.resize(dwLength);
	CryptoDestroyCipher(hRSA);
	if (rsaBlock.size() != 512) {
		fprintf(stdout,"error: RSA encryption failed\n", g_outputFile);
		exit(RET_WRITE_ERROR);
	}

	// prepare data for envelope encryption
	encryptedLicense.resize(sizeof(DWORD) + rsaBlock.size() + sizeof(DWORD) + licenseBlock.size() - (512-16));
	pData = encryptedLicense;
	*(DWORD*)pData = rsaBlock.size();
	pData += sizeof(DWORD);
	memmove(pData, rsaBlock, rsaBlock.size());
	pData += rsaBlock.size();
	*(DWORD*)pData = licenseBlock.size() - (512-16);
	pData += sizeof(DWORD);
	memmove(pData, licenseBlock + 512-16, licenseBlock.size() - (512-16));
	pData += licenseBlock.size();

	// encrypt license (envelope encryption)
	dwLength = encryptedLicense.size();
	licenseBlock.resize(16 + (dwLength & ~15) + 16);
	CryptoGetRandom(licenseBlock, 16);
	memmove(licenseBlock + 16, encryptedLicense, dwLength);
	pad = (dwLength & ~15) + 16 - dwLength;
	memset(licenseBlock + 16 + dwLength, pad, pad); 
	CryptoCreateCipher(CRYPTO_ALG_AES, &hAES);
	dwValue = CRYPTO_MODE_CBC;
	CryptoSetCipherParam(hAES, CRYPTO_CP_MODE, (BYTE*)&dwValue, sizeof(DWORD));
	Security::PRF_Expand((const BYTE*)g_HASPData, sizeof(g_HASPData), 
		(const BYTE*)&g_dwHASPID, sizeof(g_dwHASPID),
		(const BYTE*)g_HASPCrypt, sizeof(g_HASPCrypt), 
		AES_KEY, sizeof(AES_KEY));
	CryptoSetCipherParam(hAES, CRYPTO_CP_KEY, (BYTE*)AES_KEY, 16);
	CryptoSetCipherParam(hAES, CRYPTO_CP_IV, (BYTE*)AES_KEY+16, 16);
	CryptoEncrypt(hAES, NULL, licenseBlock, licenseBlock.size());
	CryptoDestroyCipher(hAES);

	// convert encrypted license to BASE64
	pDataBase64 = new char[licenseBlock.size()*2];
	dwDataBase64Length = licenseBlock.size()*2;
	CryptoConvertToBase64(licenseBlock, licenseBlock.size(), pDataBase64, &dwDataBase64Length, false);

	// write license data
	IniSetValue(hLicFile, "LicenseData", INI_SZ, pDataBase64);
	delete pDataBase64;

	fprintf(stdout,"Saving %s - ", g_outputFile);
	status = IniWriteFile(hLicFile, g_outputFile, INI_CHECKSUM_MD5_SHA1, SignLicense, NULL);
	if (status != STATUS_INI_SUCCESS) {
		fprintf(stdout,"FAILED\n", g_outputFile);
		exit(RET_WRITE_ERROR);
	}
	fprintf(stdout,"OK\n", g_outputFile);
}


bool LoadDefFile()
{
	HANDLE hDefFile;
	INISTATUS status;
	const char *pValue;
	DWORD dwValueLength;
	DWORD dwType;
	int i;

	status = IniParseFile(g_inputFile, &hDefFile, NULL, NULL);
	if (INI_FAILED(status)) {
		fprintf(stderr,"error: could not load certificate definition file %s\n", g_inputFile);
		return false;
	}

	status = IniQueryValue(hDefFile, "Application", NULL, &pValue, NULL, NULL);
	if (status != STATUS_INI_SUCCESS) {
		fprintf(stderr,"error: missing Application in %s\n", g_inputFile);
		return false;
	}
	if (strcmp(pValue, "TomoCon 3.0 Workstation") == 0) {
		g_application_type = TomoConWorkstation;
	} else if (strcmp(pValue, "TomoCon PACS Server") == 0) {
		g_application_type = PACSServer;
	} else if (strcmp(pValue, "TomoCon miniPACS Server") == 0) {
		g_application_type = miniPACSServer;
	} else {		
		fprintf(stderr,"error: unsupported application '%s' in %s\n", pValue, g_inputFile);
		return false;
	}
	
	status = IniQueryValue(hDefFile, "Type", NULL, &pValue, NULL, NULL);
	if (status != STATUS_INI_SUCCESS) {
		fprintf(stderr,"error: missing Type in %s\n", g_inputFile);
		return false;
	}
	if (strcmp(pValue, "HASP") != 0) {		
		fprintf(stderr,"error: unsupported license type '%s' in %s\n", pValue, g_inputFile);
		return false;
	}
	
	status = IniQueryValue(hDefFile, "HASPID", NULL, &pValue, NULL, NULL);
	if (status != STATUS_INI_SUCCESS) {
		fprintf(stderr,"error: missing HASPID in %s\n", g_inputFile);
		return false;
	}
	g_HASPID = strdup(pValue);
	if (sscanf(g_HASPID, "%x", &g_dwHASPID) != 1) {
		fprintf(stderr,"error: invalid format of HASPID in %s\n", g_inputFile);
		return false;
	}
	
	status = IniQueryValue(hDefFile, "HASPData", NULL, &pValue, NULL, NULL);
	if (status != STATUS_INI_SUCCESS) {
		fprintf(stderr,"error: missing HASPData in %s\n", g_inputFile);
		return false;
	}
	int d[24];
	if (24 != sscanf(pValue, "%x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x %x",
					 &d[0], &d[1], &d[2], &d[3], &d[4], &d[5], &d[6], &d[7], &d[8], &d[9], &d[10], &d[11], 
					 &d[12], &d[13], &d[14], &d[15], &d[16], &d[17], &d[18], &d[19], &d[20], &d[21], &d[22], &d[23]))
	{
		fprintf(stderr,"error: invalid format of HASPData in %s\n", g_inputFile);
		return false;
	}
	for (i = 0; i < 24; i++) 
		g_HASPData[i] = d[i];
	
	status = IniQueryValue(hDefFile, "Feature", NULL, &pValue, NULL, NULL);
	if (status != STATUS_INI_SUCCESS) {
		fprintf(stderr,"error: missing Feature in %s\n", g_inputFile);
		return false;
	}
	g_Feature = NULL;
	for (i = 0; i < sizeof(featuresTable) / sizeof(featuresTable[0]); i++) {
		if (strcmp(featuresTable[i]->featureName, pValue) == 0) {
			g_Feature = featuresTable[i];
			break;
		}
	}
	if (g_Feature == NULL) {
		fprintf(stderr,"error: unsupported feature %s in %s\n", (char*)pValue, g_inputFile);
		return false;
	}

	status = IniQueryValue(hDefFile, "FeatureData", &dwType, &pValue, &dwValueLength, NULL);
	if (status == STATUS_INI_SUCCESS) {
		g_FeatureData = new BYTE[dwValueLength];
		memmove(g_FeatureData, pValue, dwValueLength);
		g_FeatureDataSize = dwValueLength;
	} else {
		g_FeatureData = NULL;
		g_FeatureDataSize = 0;
	}

	status = IniQueryValue(hDefFile, "IssuedTo", NULL, &pValue, NULL, NULL);
	if (status != STATUS_INI_SUCCESS) {
		fprintf(stderr,"error: missing IssuedTo in %s\n", g_inputFile);
		return false;
	}
	g_IssuedTo = strdup(pValue);
	
	status = IniQueryValue(hDefFile, "ValidFrom", NULL, &pValue, NULL, NULL);
	if (status == STATUS_INI_SUCCESS && strlen(pValue) != 0)
		g_ValidFrom = strdup(pValue);
	else
		g_ValidFrom = NULL;

	status = IniQueryValue(hDefFile, "ValidTo", NULL, &pValue, NULL, NULL);
	if (status == STATUS_INI_SUCCESS && strlen(pValue) != 0)
		g_ValidTo = strdup(pValue);
	else
		g_ValidTo = NULL;
	
	status = IniQueryValue(hDefFile, "Comment", NULL, &pValue, NULL, NULL);
	if (status == STATUS_INI_SUCCESS && strlen(pValue) != 0)
		g_Comment = strdup(pValue);
	
	IniReleaseFile(hDefFile);

	sprintf(g_outputFile, "%s-%s.LIC", g_Feature->featureName, g_HASPID);
	
	return true;
}


bool LoadVendorCode(const char *path)
{
	FILE *f;
	int size;
	f = fopen(path, "rb");
	if (f == NULL)
		return false;
	size = filelength(fileno(f));
	g_HASPHL_vendor_code = new unsigned char[size];
	if (fread(g_HASPHL_vendor_code, 1, size, f) != size) {
		delete g_HASPHL_vendor_code;
		g_HASPHL_vendor_code = NULL;
		fclose(f);
		return false;
	}
	fclose(f);
	return true;
}

bool CheckHASP()
{
	if (g_application_type == TomoConWorkstation) {
		
		int p1, p2, p3, p4;
		p1 = p2 = p3 = p4 = 0;
		hasp(5, 0, 0, g_HASPPW1, g_HASPPW2, &p1, &p2, &p3, &p4);
		if (p1 != 4)
			return false;
		int port = p3;
		memmove(g_HASPCrypt, g_HASPData, sizeof(g_HASPCrypt));
		p1 = 0;
		p2 = sizeof(g_HASPCrypt);
		p3 = 0;
		p4 = (int)&g_HASPCrypt;
		hasp(61, 0, port, g_HASPPW1, g_HASPPW2, &p1, &p2, &p3, &p4);
		if (p1 != 0)
			return false;
		return true;

	} else {

		hasp_handle_t handle = 0;

		if (HASP_STATUS_OK != hasp_login(HASP_PROGNUM_DEFAULT_FID, g_HASPHL_vendor_code, &handle))
			return false;

		memmove(g_HASPCrypt, g_HASPData, sizeof(g_HASPCrypt));
		if (HASP_STATUS_OK != hasp_encrypt(handle, g_HASPCrypt, sizeof(g_HASPCrypt))) {
			hasp_logout(handle);
			return false;
		}

		hasp_logout(handle);
		return true;

	}
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

	CryptoInit();

	if (argc < 2) {
		usage();
		return RET_OK;
	}
	g_inputFile = argv[1];

	fprintf(stdout,"Loading %s\n", g_inputFile);
	if (!LoadDefFile())
		return RET_INVALID_DEF_FILE;

	if (g_application_type == TomoConWorkstation) {

		fprintf(stdout,"HASP PW1: ");
		fgets(buf, sizeof(buf), stdin);
		g_HASPPW1 = atoi(buf);

		fprintf(stdout,"HASP PW2: ");
		fgets(buf, sizeof(buf), stdin);
		g_HASPPW2 = atoi(buf);

		if (!CheckHASP()) {
			fprintf(stderr,"error: HASP not found");
			return RET_INVALID_DEF_FILE;
		}

		fprintf(stdout,"Loading tomocon-lic-rsa private key\n");
		g_tomocon_lic_rsa = ReadRSAKey("tomocon-lic-rsa.pri", false);
		if (g_tomocon_lic_rsa == NULL)
			return RET_INVALID_PEM_FILE;

		fprintf(stdout,"Loading tomocon-sig-ini-dsa private key\n");
		g_tomocon_sig_ini_dsa = ReadDSAKey("tomocon-sig-ini-dsa.pri", false);
		if (g_tomocon_sig_ini_dsa == NULL)
			return RET_INVALID_PEM_FILE;

		fprintf(stdout,"Loading tomocon-sig-lic-dsa private key\n");
		g_tomocon_sig_lic_dsa = ReadDSAKey("tomocon-sig-lic-dsa.pri", false);
		if (g_tomocon_sig_lic_dsa == NULL)
			return RET_INVALID_PEM_FILE;

	} else if (g_application_type == PACSServer || g_application_type == miniPACSServer) {

		fprintf(stdout,"HASP HL vendor code file: ");
		fgets(buf, sizeof(buf), stdin);
		char *p = buf[0] ? strchr(buf, '\n') : NULL;
		if (p)
			*p = 0;
		if (!LoadVendorCode(buf)) {
			fprintf(stderr,"error: Could not load HASP HL vendor code");
			return RET_INVALID_DEF_FILE;
		}

		if (!CheckHASP()) {
			fprintf(stderr,"error: HASP not found");
			return RET_INVALID_DEF_FILE;
		}

		fprintf(stdout,"Loading tm-pacs-lic-rsa private key\n");
		g_pacs_lic_rsa = ReadRSAKey("tm-pacs-lic-rsa.pri", false);
		if (g_pacs_lic_rsa == NULL)
			return RET_INVALID_PEM_FILE;

		fprintf(stdout,"Loading tm-pacs-sig-ini-dsa private key\n");
		g_pacs_sig_ini_dsa = ReadDSAKey("tm-pacs-sig-ini-dsa.pri", false);
		if (g_pacs_sig_ini_dsa == NULL)
			return RET_INVALID_PEM_FILE;

		fprintf(stdout,"Loading tm-pacs-sig-lic-dsa private key\n");
		g_pacs_sig_lic_dsa = ReadDSAKey("tm-pacs-sig-lic-dsa.pri", false);
		if (g_pacs_sig_lic_dsa == NULL)
			return RET_INVALID_PEM_FILE;

	}

	CreateCertificate();

	return RET_OK;	
}
