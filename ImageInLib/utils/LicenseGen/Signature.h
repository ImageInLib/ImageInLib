// Signature.h: interface for the CSignature class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_SIGNATURE_H__E221DDD9_A6A2_4633_B00D_54B3AB725A48__INCLUDED_)
#define AFX_SIGNATURE_H__E221DDD9_A6A2_4633_B00D_54B3AB725A48__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Security.h"

SECURITY_NAMESPACE_BEGIN

class CSignature  
{
public:
	CSignature();
	virtual ~CSignature();

public:
	void SetAlgo(DWORD algo);
	DWORD GetAlgo() const;

	void SetDigestType(DWORD digest);
	DWORD GetDigestType() const;

	void SetKeyID(CBinary &id);
	const CBinary& GetKeyID() const;

	void SetSignature(CBinary &sig);
	const CBinary& GetSignature() const;

protected:
	DWORD	m_Algo;
	DWORD	m_DigestType;
	CBinary	m_KeyID;
	CBinary	m_SignatureData;
};


class CSignaturesContainer
{
public:
	CSignaturesContainer();
	virtual ~CSignaturesContainer();

public:
	enum SignaturesError_T {
		errNone = 0,
		errSignatureAlreadyExists,
		errSignatureNotFound
	};

	int GetSignaturesCount() const;

	SignaturesError_T Insert(CSignature *pSignature);
	SignaturesError_T Remove(CSignature *pSignature);

	SignaturesError_T GetSignatureByKeyID(const CBinary &id, CSignature **ppSignature) const;
	SignaturesError_T GetSignatureByIndex(int idx, CSignature **ppSignature) const;

protected:
	vector<CSignature *> m_Signatures;
	map<CBinary, CSignature *> m_SignatureMap;
};

SECURITY_NAMESPACE_END

#endif // !defined(AFX_SIGNATURE_H__E221DDD9_A6A2_4633_B00D_54B3AB725A48__INCLUDED_)
