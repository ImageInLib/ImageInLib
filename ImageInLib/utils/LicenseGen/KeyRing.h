// KeyRing.h: interface for the CKeyRing class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_KEYRING_H__95F2823F_D6D7_47CD_80EA_4D6E9F83E83E__INCLUDED_)
#define AFX_KEYRING_H__95F2823F_D6D7_47CD_80EA_4D6E9F83E83E__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Security.h"

SECURITY_NAMESPACE_BEGIN

class CKeyRing  
{
public:
	CKeyRing();
	virtual ~CKeyRing();

public:
	enum KeyRingError_T {
		errNone = 0,
		errKeyAlreadyExists,
		errKeyNotFound
	};

	int GetKeysCount() const;

	KeyRingError_T Insert(CKey *pKey);
	KeyRingError_T Remove(CKey *pKey);

	KeyRingError_T GetKeyByID(const CBinary &id, CKey **ppKey) const;
	KeyRingError_T GetKeyByIndex(int idx, CKey **ppKey) const;

protected:
	vector<CKey *> m_Keys;
	map<CBinary, CKey *> m_KeyMap;
};

SECURITY_NAMESPACE_END

#endif // !defined(AFX_KEYRING_H__95F2823F_D6D7_47CD_80EA_4D6E9F83E83E__INCLUDED_)
