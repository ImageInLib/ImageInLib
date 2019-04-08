#if !defined(AFX_PRF_H__524CBE3B_3162_4F5D_AB91_A77E9DED015E__INCLUDED_)
#define AFX_PRF_H__524CBE3B_3162_4F5D_AB91_A77E9DED015E__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Security.h"

SECURITY_NAMESPACE_BEGIN

extern bool PRF_Expand(const BYTE *pSecret, DWORD cbSecret,
					   const BYTE *pLabel, DWORD cbLabel,
					   const BYTE *pSeed, DWORD cbSeed,
					   BYTE *pExpandedData, DWORD cbExpandedData);

SECURITY_NAMESPACE_END

#endif // !defined(AFX_PRF_H__524CBE3B_3162_4F5D_AB91_A77E9DED015E__INCLUDED_)
