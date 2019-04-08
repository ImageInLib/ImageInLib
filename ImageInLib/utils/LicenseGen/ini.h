/*++

Module Name:

    ini.h

Abstract:

    This include file defines the .INI file handling functions.

Author:

    Milos Opalek            (milos)

Revision History:

    Who         When        What
    --------    --------    ----------------------------------------------
    milos       21-10-98    Created.
	milos		19-12-01	API changed to return names and values by reference.
	milos		21-12-01	Ported to Linux.

Notes:

--*/

#ifndef INI_H_INCLUDED
#define INI_H_INCLUDED

typedef unsigned int INISTATUS;

#define STATUS_INI_SUCCESS                  ((INISTATUS)0)
#define STATUS_INI_NO_MORE_ENTRIES          ((INISTATUS)0x00000001)
#define STATUS_INI_CHECKSUM_OK              ((INISTATUS)0x00000002)
#define STATUS_INI_SIGNATURE_OK             ((INISTATUS)0x00000003)

#define STATUS_INI_FILE_NOT_FOUND			((INISTATUS)0x80000001)
#define STATUS_INI_FILE_READ_FAILED         ((INISTATUS)0x80000002)
#define STATUS_INI_INSUFFICIENT_RESOURCES   ((INISTATUS)0x80000003)
#define STATUS_INI_KEY_NOT_FOUND			((INISTATUS)0x80000004)
#define STATUS_INI_FILE_IS_TOO_BIG          ((INISTATUS)0x80000005)
#define STATUS_INI_INVALID_HANDLE           ((INISTATUS)0x80000006)
#define STATUS_INI_INVALID_VALUE_TYPE       ((INISTATUS)0x80000007)
#define STATUS_INI_FILE_WRITE_FAILED        ((INISTATUS)0x80000008)
#define STATUS_INI_UNSUPPORTED_CHECKSUM     ((INISTATUS)0x80000009)
#define STATUS_INI_CHECKSUM_MISSMATCH       ((INISTATUS)0x8000000A)
#define STATUS_INI_UNSUPPORTED_SIGNATURE    ((INISTATUS)0x8000000B)
#define STATUS_INI_SIGNATURE_MISSMATCH      ((INISTATUS)0x8000000C)

#define STATUS_INI_FILE_SYNTAX_ERROR        ((INISTATUS)0x90000000)
#define STATUS_INI_FILE_NOT_TOP_LEVEL       ((INISTATUS)0xa0000000)
#define STATUS_INI_FILE_KEY_EXISTS			((INISTATUS)0xb0000000)
#define STATUS_INI_FILE_EXTRA_CBRACE		((INISTATUS)0xc0000000)
#define STATUS_INI_FILE_NEWLINE				((INISTATUS)0xd0000000)
#define STATUS_INI_FILE_UNEXPECTED_EOF      ((INISTATUS)0xe0000000)

#define INI_FAILED(status)					(((status) & 0x80000000) != 0)
#define INI_SYNTAX_ERROR(status)			(((status) & 0x70000000) != 0)
#define INI_ERROR_LINE_NUMBER(status)		((status) & 0x0fffffff)

#define INI_NONE		0
#define INI_SZ			1
#define INI_MULTI_SZ	2

#define INI_CHECKSUM_NONE		0
#define INI_CHECKSUM_MD5		1
#define INI_CHECKSUM_SHA1		2
#define INI_CHECKSUM_MD5_SHA1	3

#define INI_SIGNATURE_NONE		0
#define INI_SIGNATURE_RSA		1
#define INI_SIGNATURE_DSA		2

typedef INISTATUS (*INI_MAKE_SIGNATURE_CALLBACK)(
	VOID       *pCallbackCtx,
	DWORD        dwDigestType,
	const BYTE  *pDigest,
	DWORD		 cbDigest,
	DWORD       *pdwSigType,
	BYTE       **ppPubKeyID,	// must be allocated by malloc !
	DWORD		*pcbPubKeyID,
	BYTE       **ppSignature,   // must be allocated by malloc !
	DWORD		*pcbSignature
	);

typedef INISTATUS (*INI_VERIFY_SIGNATURE_CALLBACK)(
	VOID       *pCallbackCtx,
	DWORD       dwDigestType,
	const BYTE *pDigest,
	DWORD		cbDigest,
	DWORD       dwSigType,
	const BYTE *pPubKeyID,
	DWORD		cbPubKeyID,
	const BYTE *pSignature,
	DWORD		cbSignature
	);

extern INISTATUS
IniParseFile(
    const char *lpFileName,
    HANDLE *hIniFile,
	INI_VERIFY_SIGNATURE_CALLBACK pVerifyCallback = NULL,
	void *pVerifyCallbackCtx = NULL
    );

extern INISTATUS
IniParseBuffer(
    const char *lpBuffer,
    HANDLE *hIniFile,
	INI_VERIFY_SIGNATURE_CALLBACK pVerifyCallback = NULL,
	void *pVerifyCallbackCtx = NULL
    );

extern VOID
IniReleaseFile(
    HANDLE hIniFile
    );

extern INISTATUS
IniOpenKey(
    HANDLE       hKey,
    const char  *lpSubKey,
    HANDLE      *phResult
    );

extern INISTATUS
IniEnumKey(
    HANDLE       hKey,
    DWORD        dwIndex,
    const char **lpName,
	DWORD       *pLineNumber = NULL
    );
 
extern INISTATUS
IniEnumValue(
    HANDLE       hKey,
    DWORD        dwIndex,
    const char **ppName,
    DWORD       *pdwType,
    const char **pData,
    DWORD       *pcbData,
	DWORD       *pLineNumber = NULL
    );
 
extern INISTATUS
IniQueryInfoKey(
    HANDLE  hKey,
    DWORD  *pcSubKeys,
    DWORD  *pcValues,
	DWORD  *pLineNumber = NULL
    );

extern INISTATUS
IniQueryValue(
    HANDLE       hKey,
    const char  *pValueName,
    DWORD       *pdwType,
    const char **pData,
    DWORD       *pcbData,
	DWORD       *pLineNumber = NULL
    );


extern INISTATUS 
IniCreateFile(
	HANDLE* phIniFile
	);

extern INISTATUS
IniWriteFile(
    HANDLE hIniFile,
    const char *pFileName,
	DWORD dwChecksumType = INI_CHECKSUM_NONE,
	INI_MAKE_SIGNATURE_CALLBACK pSignCallback = NULL,
	void *pSignCallbackCtx = NULL
    );

extern INISTATUS
IniWriteBuffer(
    HANDLE hIniFile,
    string &strBuffer,
	DWORD dwChecksumType = INI_CHECKSUM_NONE,
	INI_MAKE_SIGNATURE_CALLBACK pSignCallback = NULL,
	void *pSignCallbackCtx = NULL
    );

extern INISTATUS IniCreateKey(
	HANDLE		hKey,
	const char *pKeyName,
	HANDLE	   *phkResult
	);

extern INISTATUS
IniSetValue(
	HANDLE	    hKey,
	const char *pValueName,
	DWORD       dwType,
	const char *lpData
	);

extern INISTATUS
IniDeleteKey(
    HANDLE hKey,
	const char *pSubName
    );

extern INISTATUS
IniDeleteValue(
    HANDLE hKey,
	const char *pValueName
    );

extern const char *
IniErrorMessage(
	INISTATUS status
	);
 
#endif /* INI_H_INCLUDED */
