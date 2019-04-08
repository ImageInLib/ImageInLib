// PRF.cpp: implementation of the PRF_Expand.
//
//////////////////////////////////////////////////////////////////////

#include "headers.h"
#include "PRF.h"

using namespace ::Security;

bool Security::PRF_Expand(
	const BYTE *pSecret, DWORD cbSecret,
	const BYTE *pLabel, DWORD cbLabel,
	const BYTE *pSeed, DWORD cbSeed,
	BYTE *pExpandedData, DWORD cbExpandedData)
{
	HCRYPTOHASH hHMAC = NULL;
	const BYTE *pS1;
	DWORD L_S1;
	const BYTE *pS2;
	DWORD L_S2;
	BYTE HMAC[20];
	BYTE A[20];
	int pass;
	BYTE *pBlock;

	if (cbSecret == 0)
		goto fail;

	// get S1 and S2
	L_S1 = L_S2 = (cbSecret + 1) / 2;
	pS1 = pSecret;
	pS2 = pSecret + cbSecret - L_S2;

	for (pass = 0; pass < 2; pass++) {
		DWORD dwHMACLen = pass ? 20 : 16;
		DWORD dwLen;
		DWORD dwBlock;

		// initialize HMAC
		if (!CryptoCreateHash(pass ? CRYPTO_ALG_HMAC_SHA1 : CRYPTO_ALG_HMAC_MD5, &hHMAC))
			goto fail;
		if (!CryptoSetHashParam(hHMAC, CRYPTO_HP_KEY, pass ? pS2 : pS1, pass ? L_S2 : L_S2))
			goto fail;

		// calculate A(1)
		if (!CryptoHashData(hHMAC, pLabel, cbLabel))
			goto fail;
		if (!CryptoHashData(hHMAC, pSeed, cbSeed))
			goto fail;
		dwLen = dwHMACLen;
		if (!CryptoGetHashParam(hHMAC, CRYPTO_HP_HASHVALUE, A, &dwLen))
			goto fail;

		for (dwBlock = 0; dwBlock < (cbExpandedData + dwHMACLen - 1) / dwHMACLen; dwBlock++) {

			// HMAC_hash(secret, A(i) + seed)
			if (!CryptoHashData(hHMAC, A, dwHMACLen))
				goto fail;
			if (!CryptoHashData(hHMAC, pLabel, cbLabel))
				goto fail;
			if (!CryptoHashData(hHMAC, pSeed, cbSeed))
				goto fail;
			dwLen = dwHMACLen; 
			if (!CryptoGetHashParam(hHMAC, CRYPTO_HP_HASHVALUE, HMAC, &dwLen))
				goto fail;

			dwLen = min(dwHMACLen, cbExpandedData - dwHMACLen * dwBlock);
			pBlock = pExpandedData + dwHMACLen * dwBlock;
			if (pass == 0) {
				memmove(pBlock, HMAC, dwLen);
			} else {
				for (int i = 0; i < (int)dwLen; i++)
					pBlock[i] ^= HMAC[i];
			}

			// A(i) = HMAC_hash(secret, A(i-1))
			if (!CryptoHashData(hHMAC, A, dwHMACLen))
				goto fail;
			dwLen = dwHMACLen;
			if (!CryptoGetHashParam(hHMAC, CRYPTO_HP_HASHVALUE, A, &dwLen))
				goto fail;
		}

		// destroy HMAC
		CryptoDestroyHash(hHMAC);
		hHMAC = NULL;
	}

	memset(A, 0, sizeof(A));
	memset(HMAC, 0, sizeof(HMAC));
	return true;

fail:
	memset(A, 0, sizeof(A));
	memset(HMAC, 0, sizeof(HMAC));
	if (hHMAC)
		CryptoDestroyHash(hHMAC);
	memset(pExpandedData, 0, cbExpandedData);
	return false;
}
