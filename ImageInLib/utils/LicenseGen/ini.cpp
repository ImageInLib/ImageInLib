/*++

Module Name:

    ini.cpp

Abstract:

    This module contains the .INI files handling code.

Author:

    Milos Opalek            (milos)

Revision History:

    Who			When		What
    --------	--------	---------------------------------------------------
    milos		21-10-98    Created.
	milos		19-12-01	API changed to return names and values by reference.
	milos		21-12-01	Ported to Linux.
	milos		08-08-04	Digital signatures support.

Notes:
	thread safe.

--*/

#include "headers.h"
#include "../crypto/crypto.h"
#include "ini.h"
#include <fcntl.h>
#include <sys/stat.h>
#include <stdio.h>

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif


#define CHECKSUM_BEGIN			"\n#----CHECKSUM BEGIN----\n"
#define CHECKSUM_END			"\n#-----CHECKSUM END-----\n"
#define CHECKSUM_HASH_MD5		"# Hash:MD5\n#\n# "
#define CHECKSUM_HASH_SHA1		"# Hash:SHA1\n#\n# "
#define CHECKSUM_HASH_MD5_SHA1	"# Hash:MD5+SHA1\n#\n# "

#define SIGNATURE_BEGIN			"\n#----SIGNATURE BEGIN----\n"
#define SIGNATURE_END			"\n#-----SIGNATURE END-----\n"
#define SIGNATURE_RSA			"# Signature:RSA\n#\n# "
#define SIGNATURE_DSA			"# Signature:DSA\n#\n# "


//
// Type definitions.
//

typedef struct ININODE {
    DWORD    Magic;
    DWORD    Type;
    char    *Name;
    DWORD    ValueSize;
	DWORD    LineNumber;
	bool	 allocatedValue;
    union {
        char           *Value;
        struct ININODE *SubKey;
    };
    ININODE *Parent; 
    ININODE *Next; 
} ININODE, *PININODE;

typedef struct INIFILE {
    DWORD    Magic;
    char*    Buffer;
    DWORD    BufferSize;
    DWORD    BufferPos;
    ININODE  RootNode;
} INIFILE, *PINIFILE;

#define INIFILE_MAGIC ('I'+('N'<<8)+('I'<<16)+('F'<<24))
#define ININODE_MAGIC ('I'+('N'<<8)+('I'<<16)+('N'<<24))
#define INISECT_MAGIC ('I'+('N'<<8)+('I'<<16)+('S'<<24))

typedef struct BUF_STATE {
	string *stringBuffer;
	int  fileHandle;
	HCRYPTOHASH hMD5;
	HCRYPTOHASH hSHA1;
} BUF_STATE;


//
// Function prototypes.
//

static char * 
IniReadLine(
    PINIFILE IniFile
    );

static char *
IniParseQuotedString(
    char *String
    );

static INISTATUS
IniWriteNode(
	BUF_STATE *pBufState,
	ININODE *pNode,
	int indentLevel);

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


char *
IniReadLine(
    PINIFILE IniFile
    )

/*++

Routine Description:

    This routine reads one line from the .INI file and strips comments and
    leading and trailing whitespaces.

Arguments:

Return Value:

--*/

{
    char *p, *lineStart;
    DWORD pos;
    DWORD state;

    pos = IniFile->BufferPos;
    p = lineStart = IniFile->Buffer + IniFile->BufferPos;
    state = 0;
    
    if (pos >= IniFile->BufferSize) {
        return NULL;
    }

    for (;;) {
        char c;

        c = IniFile->Buffer[pos];

        if (c == '\n' || pos == IniFile->BufferSize) {
            do --p; while (p > IniFile->Buffer && (*p == '\t' || *p == ' '));
            *++p = 0;
            IniFile->BufferPos = pos + 1;
            return lineStart;
        }

        if (c == '\r') {
            pos++;
            continue;
        }

        switch (state) {

        case 0: // skip leading white spaces
            if (c == '\t' || c == ' ')
				break;
            state = 1;
        case 1: // normal string
            if (c == '\\') {
                state = 4;
            } else if (c == '"') {
                state = 3;
            } else if (c == '#') {
                state = 2;
                break;
            }
            *p++ = c;
            break;

        case 2: // comment
            break;

        case 3: // quoted string
            if (c == '\\') {
                state = 4;
            } else if (c == '"') {
                state = 0;
            }
            *p++ = c;
            break;

        case 4: // escape sequences
            *p++ = c;
            state = 3;
            break;

        }

        pos++;
    }

} // IniReadLine



char *
IniParseQuotedString(
    char *String
    )

/*++

Routine Description:

    This routine parses quoted string.

Arguments:

Return Value:

--*/

{
    char *p, *s, *e;
    DWORD state;
    DWORD num = 0;
    char c;

    p = s = e = String;
    state = 0;
    
    for (;;) {

        c = *s;

        if (c == '\0') {
            *e = '\0';
            return String;
        }

        switch (state) {

        case 0: // skip leading white spaces
            if (c != '\t' && c != ' ') {
				if (c == '"') {
					state = 2;
					break;
				}
                state = 1;
                *p++ = c;
                e = p;
            }
            break;

        case 1: // normal string
            if (c == '"') {
                state = 2;
                break;
            }
            *p++ = c;
            if (c != '\t' && c != ' ') {
                e = p;
            }
            break;

        case 2: // quoted string
            if (c == '\\') {
                state = 3;
                break;
            } else if (c == '"') {
                state = 1;
                break;
            }
            *p++ = c;
            e = p;
            break;

        case 3: // escape sequences
            if (c == 'a') {
                c = '\a';
            } else if (c == 'b') {
                c = '\b';
            } else if (c == 'f') {
                c = '\f';
            } else if (c == 'n') {
                c = '\n';
            } else if (c == 'r') {
                c = '\r';
            } else if (c == 't') {
                c = '\t';
            } else if (c == 'v') {
                c = '\v';
            } else if (c == '"') {
                c = '"';
            } else if (c == '\\') {
                c = '\\';
            } else if (c >= '0' && c <= '7') {
                num = c - '0';
                state = 4;
                break;
            } else  if (c == 'x') {
				num = 0;
                state = 5;
                break;
            }
            *p++ = c;
			e = p;
            state = 2;
            break;

        case 4: // ASCII character in octal notation
            if (c >= '0' && c <= '7') {
                num = num * 8 + (c - '0');
                break;
            }
            *p++ = (char)num;
			num = 0;
			e = p;
            state = 2;
            break;

        case 5: // ASCII character in hexadecimal notation
            if ((c >= '0' && c <= '9') || 
                (c >= 'A' && c <= 'F') ||
                (c >= 'a' && c <= 'f')) 
            {
                if (c >= '0' && c <= '9') {
                    c -= '0';
                } else {
                    c = (c | 0x20) - 'a';
                }
                num = num * 16 + c;
                if (num < 0x10) {
                    break;
                }
            }
            *p++ = (char)num;
			num = 0;
			e = p;
            state = 2;
            break;
        }

        s++;
    }

} // ParseQuotedString



INISTATUS
IniParseFile(
    const char *pFileName,
    HANDLE *hIniFile,
	INI_VERIFY_SIGNATURE_CALLBACK pVerifyCallback,
	void *pVerifyCallBackCtx
    )

{
    INISTATUS status;
    int handle = -1;
    char *buffer = NULL;
    DWORD iniFileSize;
	
    //
    // Open .INI file.
    //
#ifdef _WIN32
	handle = open(pFileName, O_RDONLY | O_BINARY);
#else
	handle = open(pFileName, O_RDONLY);
#endif

    if (handle == -1) {
        status = STATUS_INI_FILE_NOT_FOUND;
        goto fail;
    }

    //
    // Allocate buffer.
    //
    iniFileSize = lseek(handle, 0, SEEK_END);
	lseek(handle, 0, SEEK_SET);
    buffer = (char*)malloc(iniFileSize + 1);
    if (buffer == NULL) {
        status = STATUS_INI_INSUFFICIENT_RESOURCES;
        goto fail;
    }

    //
    // Read file.
    //
    int bytesRead;
    if ((bytesRead = read(handle, buffer, iniFileSize)) == -1) {
        status = STATUS_INI_FILE_READ_FAILED;
        goto fail;
    }
	buffer[bytesRead] = 0;

	//
	// Close file.
	//
    close(handle);
	handle = -1;

	status = IniParseBuffer(buffer, hIniFile, pVerifyCallback, pVerifyCallBackCtx);
	free(buffer);

    return status;

fail:
	if (handle != -1)
		close(handle);
	if (buffer)
		free(buffer);

    return status;

}


INISTATUS
IniParseBuffer(
    const char *lpBuffer,
    HANDLE *hIniFile,
	INI_VERIFY_SIGNATURE_CALLBACK pVerifyCallback,
	void *pVerifyCallbackCtx
    )

{
    INISTATUS status;
    PINIFILE iniFile = NULL;
    char *buffer = NULL;
    char *line, *p;
    DWORD keyLine = 0;
    DWORD lineNumber;
    DWORD parseState;
    PININODE parentNode;
    PININODE currentNode;
    DWORD iniFileSize;
	DWORD textLen;
	HCRYPTOHASH hMD5 = NULL;
	HCRYPTOHASH hSHA1 = NULL;
	DWORD dwHashType = INI_CHECKSUM_NONE;
	DWORD dwSigType = INI_SIGNATURE_NONE;
	BYTE *pPubKeyID = NULL;
	DWORD cbPubKeyID;
	BYTE *pSignature = NULL;
	DWORD cbSignature;
	char *sig = 0;
	char *chksum = 0;
	bool bChecksum = false;
	bool bSignature = false;
	int nHashLen = 16+20;
	char hash[16+20];

    //
    // Copy buffer.
    //
    iniFileSize = strlen(lpBuffer);
    buffer = (char*)malloc(iniFileSize + 1);
    if (buffer == NULL) {
        status = STATUS_INI_INSUFFICIENT_RESOURCES;
        goto fail;
    }
	strcpy(buffer, lpBuffer);

	//
	// Check file checksum
	//
	chksum = strstr(buffer, CHECKSUM_BEGIN);
	while (chksum) {
		p = strstr(chksum + 1, CHECKSUM_BEGIN);
		if (p == NULL)
			break;
		chksum = p;
	}
	if (chksum != NULL) {
		char *pHashIDStr;
		char hashhex[256];
		p = chksum + strlen(CHECKSUM_BEGIN);
		if (strncmp(p, CHECKSUM_HASH_MD5, strlen(CHECKSUM_HASH_MD5)) == 0) {
			dwHashType = INI_CHECKSUM_MD5;
			pHashIDStr = CHECKSUM_HASH_MD5;
			nHashLen = 16;
			if (!CryptoCreateHash(CRYPTO_ALG_MD5, &hMD5)) {
				status = STATUS_INI_INSUFFICIENT_RESOURCES;
				goto fail;
			}
		} else if (strncmp(p, CHECKSUM_HASH_SHA1, strlen(CHECKSUM_HASH_SHA1)) == 0) {
			dwHashType = INI_CHECKSUM_SHA1;
			pHashIDStr = CHECKSUM_HASH_SHA1;
			nHashLen = 20;
			if (!CryptoCreateHash(CRYPTO_ALG_SHA1, &hSHA1)) {
				status = STATUS_INI_INSUFFICIENT_RESOURCES;
				goto fail;
			}
		} else if (strncmp(p, CHECKSUM_HASH_MD5_SHA1, strlen(CHECKSUM_HASH_SHA1)) == 0) {
			dwHashType = INI_CHECKSUM_MD5_SHA1;
			pHashIDStr = CHECKSUM_HASH_MD5_SHA1;
			nHashLen = 16+20;
			if (!CryptoCreateHash(CRYPTO_ALG_MD5, &hMD5)) {
				status = STATUS_INI_INSUFFICIENT_RESOURCES;
				goto fail;
			}
			if (!CryptoCreateHash(CRYPTO_ALG_SHA1, &hSHA1)) {
				status = STATUS_INI_INSUFFICIENT_RESOURCES;
				goto fail;
			}
		} else {
			status = STATUS_INI_UNSUPPORTED_CHECKSUM;
			goto fail;
		}
		p += strlen(pHashIDStr);
		textLen = chksum - buffer ;
		if (hSHA1 && !CryptoHashData(hSHA1, (BYTE*)buffer, textLen)) {
			status = STATUS_INI_INSUFFICIENT_RESOURCES;
			goto fail;
		}
		if (hMD5 && !CryptoHashData(hMD5, (BYTE*)buffer, textLen)) {
			status = STATUS_INI_INSUFFICIENT_RESOURCES;
			goto fail;
		}
		textLen = 16;
		if (hMD5 && !CryptoGetHashParam(hMD5, CRYPTO_HP_HASHVALUE, (BYTE*)hash, &textLen)) {
			status = STATUS_INI_INSUFFICIENT_RESOURCES;
			goto fail;
		}
		textLen = 20;
		if (hSHA1 && !CryptoGetHashParam(hSHA1, CRYPTO_HP_HASHVALUE, (BYTE*)hash + (hMD5 ? 16 : 0), &textLen)) {
			status = STATUS_INI_INSUFFICIENT_RESOURCES;
			goto fail;
		}
		if (hMD5)
			CryptoDestroyHash(hMD5);
		hMD5 = NULL;
		if (hSHA1)
			CryptoDestroyHash(hSHA1);
		hSHA1 = NULL;
		textLen = sizeof(hashhex);
		if (!CryptoConvertToHex((BYTE*)hash, nHashLen, hashhex, &textLen)) {
			status = STATUS_INI_INSUFFICIENT_RESOURCES;
			goto fail;
		}
		strcat(hashhex, CHECKSUM_END);
		textLen = strlen(hashhex);
		if (strncmp(hashhex, p, textLen) != 0) {
			status = STATUS_INI_CHECKSUM_MISSMATCH;
			goto fail;
		}
		if (strncmp(p + textLen, SIGNATURE_BEGIN, strlen(SIGNATURE_BEGIN)) == 0) {
			sig = p + textLen;
		} else {
			if (p[textLen] != '\0') {
				status = STATUS_INI_CHECKSUM_MISSMATCH;
				goto fail;
			}
			sig = NULL;
		}
		bChecksum = true;
	} else {
		bChecksum = false;
	}

	//
	// Get file signature
	//
	if (sig == NULL) {
		sig = strstr(buffer, SIGNATURE_BEGIN);
		while (sig) {
			p = strstr(sig + 1, SIGNATURE_BEGIN);
			if (p == NULL)
				break;
			sig = p;
		}
	}
	if (sig != NULL && pVerifyCallback) {
		char *pSigIDStr;
		char *pend;
		if (!bChecksum) {
			status = STATUS_INI_SIGNATURE_MISSMATCH;
			goto fail;
		}
		p = sig + strlen(SIGNATURE_BEGIN);
		if (strncmp(p, SIGNATURE_RSA, strlen(SIGNATURE_RSA)) == 0) {
			dwSigType = INI_SIGNATURE_RSA;
			pSigIDStr = SIGNATURE_RSA;
		} else if (strncmp(p, SIGNATURE_DSA, strlen(SIGNATURE_DSA)) == 0) {
			dwSigType = INI_SIGNATURE_RSA;
			pSigIDStr = SIGNATURE_RSA;
		} else {
			status = STATUS_INI_UNSUPPORTED_SIGNATURE;
			goto fail;
		}
		p += strlen(pSigIDStr);
		// get public key identifier
		if (*p == '[') {
			char *pend = strchr(p, ']');
			if (pend == NULL) {
				status = STATUS_INI_SIGNATURE_MISSMATCH;
				goto fail;
			}
			p++;
			*pend = '\0';
			textLen = pend - p;
			if (!CryptoConvertFromBase64(p, NULL, &cbPubKeyID, false)) {
				*pend = ']';
				status = STATUS_INI_SIGNATURE_MISSMATCH;
				goto fail;
			}
			pPubKeyID = (BYTE*)malloc(cbPubKeyID);
			if (pPubKeyID == NULL) {
				*pend = ']';
				status = STATUS_INI_SIGNATURE_MISSMATCH;
				goto fail;
			}
			if (!CryptoConvertFromBase64(p, pPubKeyID, &cbPubKeyID, false)) {
				*pend = ']';
				status = STATUS_INI_SIGNATURE_MISSMATCH;
				goto fail;
			}
			*pend = ']';
			p = pend+1;
		}
		// get signature
		pend = strchr(p, '\n');
		if (pend == NULL) {
			status = STATUS_INI_SIGNATURE_MISSMATCH;
			goto fail;
		}
		*pend = '\0';
		textLen = pend - p;
		if (!CryptoConvertFromBase64(p, NULL, &cbSignature, false)) {
			*pend = '\n';
			status = STATUS_INI_SIGNATURE_MISSMATCH;
			goto fail;
		}
		pSignature = (BYTE*)malloc(cbSignature);
		if (pSignature == NULL) {
			*pend = '\n';
			status = STATUS_INI_SIGNATURE_MISSMATCH;
			goto fail;
		}
		if (!CryptoConvertFromBase64(p, pSignature, &cbSignature, false)) {
			*pend = '\n';
			status = STATUS_INI_SIGNATURE_MISSMATCH;
			goto fail;
		}
		*pend = '\n';
		p = pend;
		if (strcmp(p, SIGNATURE_END) != 0) {
			status = STATUS_INI_SIGNATURE_MISSMATCH;
			goto fail;
		}
		// verify signature
		status = pVerifyCallback(pVerifyCallbackCtx, dwHashType, (BYTE*)&hash, 
			nHashLen, dwSigType, pPubKeyID, cbPubKeyID, pSignature, cbSignature);
		if (INI_FAILED(status))
			goto fail;
		free(pPubKeyID);
		pPubKeyID = NULL;
		free(pSignature);
		pSignature = NULL;
		bSignature = true;
	} else {
		bSignature = false;
	}

    //
    // Allocate and initialize INIFILE structure.
    //
    iniFile = (PINIFILE)malloc(sizeof(INIFILE));
    if (iniFile == NULL) {
        status = STATUS_INI_INSUFFICIENT_RESOURCES;
        goto fail;
    }
    iniFile->Magic = INIFILE_MAGIC;
    iniFile->Buffer = buffer;
    iniFile->BufferSize = iniFileSize;
    iniFile->BufferPos = 0;
    iniFile->RootNode.Magic = ININODE_MAGIC;
    iniFile->RootNode.Type = INI_NONE;
    iniFile->RootNode.Next = NULL;
    iniFile->RootNode.Parent = NULL;
    iniFile->RootNode.Name = NULL;
    iniFile->RootNode.ValueSize = 0;
    iniFile->RootNode.Value = NULL;
	iniFile->RootNode.allocatedValue = false;
	iniFile->RootNode.SubKey = NULL;

    //
    // Parse file.
    //
    lineNumber = 1;
    parseState = 0;
    currentNode = &iniFile->RootNode;
    parentNode = &iniFile->RootNode;
    parentNode->SubKey = NULL;

	do {
		line = IniReadLine(iniFile);
	} while (line && *line == '\0');
    DWORD multiLen;
    char *multi;
    DWORD nameLen;
    char *name;
	multiLen = 0; multi = NULL;
	nameLen = 0; name = NULL;
    while (line != NULL || parseState == 6) {
        PININODE newNode, node;

        switch (parseState) {

        case 0:

            if (*line == '[') {
                //
                // Section header.
                //

                if (parentNode != &iniFile->RootNode) {
                    //
                    // Section header not at top level.
                    //
                    status = STATUS_INI_FILE_NOT_TOP_LEVEL + (lineNumber & 0xffffff);
                    goto fail;
                }
                if (line[strlen(line)-1] != ']') {
                    //
                    // Missing ']' in section header.
                    //
                    status = STATUS_INI_FILE_SYNTAX_ERROR + (lineNumber & 0xffffff);
                    goto fail;
                }

                //
                // Get key name.
                //
                line[strlen(line)-1] = '\0';
                line++;
                IniParseQuotedString(line);

                //
                // Add new key
                //
                newNode = (PININODE)malloc(sizeof(ININODE) + strlen(line) + 1);
                if (newNode == NULL) {
                    status = STATUS_INI_INSUFFICIENT_RESOURCES;
                    goto fail;
                }
                newNode->Magic = INISECT_MAGIC;
                newNode->Type = INI_NONE;
				newNode->LineNumber = lineNumber;
                newNode->Name = (char*)(newNode + 1);
                newNode->Value = NULL;
                newNode->Parent = parentNode;
				newNode->allocatedValue = false;
                strcpy(newNode->Name, line);
				bool bAlreadyExists = false;
                for (node = parentNode->SubKey; node != NULL; node = node->Next) {
                    if (strcmp(node->Name, newNode->Name) == 0) {
						bAlreadyExists = true;
						break;
                    }
                    if (node->Next == NULL) {
                        break;
                    }
                }
				if (bAlreadyExists) {
					free(newNode);
					currentNode = node;
				} else {
					newNode->Next = NULL;
					if (node == NULL) {
						parentNode->SubKey = newNode;
					} else {
						node->Next = newNode;
					}
					currentNode = newNode;
				}

                *line = '\0';
                break;
            }

            if (*line == '}') {
                //
                // Closing brace.
                //
                line++; while (*line == ' ' || *line == '\t') line++;
                if (*line != '\0') {
                    status = STATUS_INI_FILE_SYNTAX_ERROR + (lineNumber & 0xffffff);
                    goto fail;
                }
                if (parentNode == &iniFile->RootNode && currentNode->Magic == INISECT_MAGIC) {
                    status = STATUS_INI_FILE_EXTRA_CBRACE + (lineNumber & 0xffffff);
                    goto fail;
                }
                currentNode = parentNode;
                parentNode = parentNode->Parent;
				if (parentNode == NULL)
					parentNode = &iniFile->RootNode;
                break;
            }

            if (*line == '=' || *line == '{') {
                status = STATUS_INI_FILE_SYNTAX_ERROR + (lineNumber & 0xffffff);
                goto fail;
            }

            //
            // Subkey or relation.
            //
            for (name = line, nameLen = 0; *line && *line != '='; nameLen++, line++) 
            {
                if (*line == '=') {
                    while (name[nameLen-1] == ' ' || name[nameLen-1] == '\t') {
                        nameLen--;
                    }
                    break;
                }
                if (*line == '"') {
                    nameLen++; line++;
                    while (*line && *line != '"') {
                        if (*line == '\\' && *(line + 1)) {
                            nameLen++; line++;
                        }
                        nameLen++; line++;
                    }
                    if (*line == '\0') {
                        //
                        // Missing closing quotation mark.
                        //
                        status = STATUS_INI_FILE_NEWLINE + (lineNumber & 0xffffff);
                        goto fail;
                    }
                }
            }

            if (*line == '\0') {
                //
                // Missing '=', multi-string
                //
                line = name;
                parseState = 5;
                continue;
            }
            line++; while (*line == ' ' || *line == '\t') line++;
            keyLine = lineNumber;
            parseState = 1;

            break;

        case 1:
            //
            // Expecting '{' or string.
            //
            if (*line == '{') {
                //
                // Subkey or multi-string.
                //
                line++; while (*line == ' ' || *line == '\t') line++;
                if (*line != '\0') {
                    status = STATUS_INI_FILE_SYNTAX_ERROR + (lineNumber & 0xffffff);
                    goto fail;
                }
                parseState = 2;
                break;
            }

            //
            // String value.
            //
            for (p = line; *p; p++) 
            {
                if (*p == '"') {
                    p++;
                    while (*p && *p != '"') {
                        if (*p == '\\' && *(p + 1)) {
                            p++;
                        }
                        p++;
                    }
                    if (*p == '\0') {
                        //
                        // Missing closing quotation mark.
                        //
                        status = STATUS_INI_FILE_NEWLINE + (lineNumber & 0xffffff);
                        goto fail;
                    }
                }
            }

			if (name) {
				name[nameLen] = '\0';
				IniParseQuotedString(name);
			}
            IniParseQuotedString(line);

            //
            // Add new string.
            //
            newNode = (PININODE)malloc(sizeof(ININODE) + strlen(name) + 1 + strlen(line) + 1);
            if (newNode == NULL) {
                status = STATUS_INI_INSUFFICIENT_RESOURCES;
                goto fail;
            }
            newNode->Magic = ININODE_MAGIC;
            newNode->Type = INI_SZ;
            newNode->Name = (char*)(newNode + 1);
			newNode->LineNumber = keyLine;
            newNode->ValueSize = strlen(line) + 1;
            newNode->Value = newNode->Name + strlen(name) + 1;
            newNode->Parent = currentNode;
			newNode->allocatedValue = false;
            strcpy(newNode->Name, name);
            strcpy((char*)newNode->Value, line);
            for (node = currentNode->SubKey; node != NULL; node = node->Next) {
                if (strcmp(node->Name, newNode->Name) == 0) {
                    status = STATUS_INI_FILE_KEY_EXISTS + (keyLine & 0xffffff);
                    goto fail;
                }
                if (node->Next == NULL) {
                    break;
                }
            }
            newNode->Next = NULL;
            if (node == NULL) {
                currentNode->SubKey = newNode;
            } else {
                node->Next = newNode;
            }

            *line = '\0';
            parseState = 0;
            break;

        case 2:
            //
            // Subkey or multi-string.
            //
            parseState = 3;
            for (p = line; *p; p++) 
            {
                if (*p == '=') {
                    //
                    // Add subkey.
                    //
                    name[nameLen] = '\0';
                    IniParseQuotedString(name);

                    newNode = (PININODE)malloc(sizeof(ININODE) + strlen(name) + 1);
                    if (newNode == NULL) {
                        status = STATUS_INI_INSUFFICIENT_RESOURCES;
                        goto fail;
                    }
                    newNode->Magic = ININODE_MAGIC;
                    newNode->Type = INI_NONE;
					newNode->LineNumber = keyLine;
                    newNode->Name = (char*)(newNode + 1);
                    newNode->Value = NULL;
                    newNode->Parent = currentNode;
					newNode->allocatedValue = false;
                    strcpy((char*)newNode->Name, name);
                    for (node = currentNode->SubKey; node != NULL; node = node->Next) {
                        if (strcmp(node->Name, newNode->Name) == 0) {
                            status = STATUS_INI_FILE_KEY_EXISTS + (keyLine & 0xffffff);
                            goto fail;
                        }
                        if (node->Next == NULL) {
                            break;
                        }
                    }
                    newNode->Next = NULL;
                    if (node == NULL) {
                        currentNode->SubKey = newNode;
                    } else {
                        node->Next = newNode;
                    }
                    parentNode = currentNode;
                    currentNode = newNode;

                    parseState = 0;
                    break;
                }

                if (*p == '"') {
                    p++;
                    while (*p && *p != '"') {
                        if (*p == '\\' && *(p + 1)) {
                            p++;
                        }
                        p++;
                    }
                    if (*p == '\0') {
                        //
                        // Missing closing quotation mark.
                        //
                        status = STATUS_INI_FILE_NEWLINE + (lineNumber & 0xffffff);
                        goto fail;
                    }
                }
            }
            break;

        case 3:
            //
            // Multi-string.
            //
            if (*line == '}' && *(line+1) == '\0') {
                multi = "";
                multiLen = 1;
                parseState = 4;
                break;
            }

            IniParseQuotedString(line);
            multiLen = strlen(line) + 1;
            multi = line;
            line = "";

            parseState = 4;
            break;

        case 4:
            //
            // Multi-string.
            //
            if (*line == '}' && *(line+1) == '\0') {
                //
                // Add multistring.
                //
                name[nameLen] = '\0';
                IniParseQuotedString(name);

                newNode = (PININODE)malloc(sizeof(ININODE) + 
                                                  strlen(name) + 1 + 
                                                  multiLen + 1);
                if (newNode == NULL) {
                    status = STATUS_INI_INSUFFICIENT_RESOURCES;
                    goto fail;
                }
                newNode->Magic = ININODE_MAGIC;
                newNode->Type = INI_MULTI_SZ;
				newNode->LineNumber = keyLine;
                newNode->Name = (char*)(newNode + 1);
                newNode->ValueSize = multiLen + 1;
                newNode->Value = newNode->Name + strlen(name) + 1;
                newNode->Parent = currentNode;
				newNode->allocatedValue = false;
                strcpy(newNode->Name, name);
                memmove(newNode->Value, multi, multiLen);
                newNode->Value[multiLen] = 0;
                for (node = currentNode->SubKey; node != NULL; node = node->Next) {
                    if (strcmp(node->Name, newNode->Name) == 0) {
                        status = STATUS_INI_FILE_KEY_EXISTS + (keyLine & 0xffffff);
                        goto fail;
                    }
                    if (node->Next == NULL) {
                        break;
                    }
                }
                newNode->Next = NULL;
                if (node == NULL) {
                    currentNode->SubKey = newNode;
                } else {
                    node->Next = newNode;
                }

                line = "";
                parseState = 0;
                break;
            }

            IniParseQuotedString(line);
            memmove(multi + multiLen, line, strlen(line)+1);
            multiLen += strlen(multi + multiLen) + 1;
            line = "";

            break;

        case 5:
            //
            // Multi-string.
            //
            if (*line == '[') {
                multi = "";
                multiLen = 1;
                parseState = 6;
                continue;
            }

            IniParseQuotedString(line);
            multiLen = strlen(line) + 1;
            multi = line;
            line = "";

            parseState = 6;
            break;

        case 6:
            //
            // Multi-string.
            //
			if (line == NULL && currentNode == &iniFile->RootNode) {
				status = STATUS_INI_FILE_UNEXPECTED_EOF + (lineNumber & 0xffffff);
				goto fail;
			}
            if (line == NULL || *line == '[') {
                //
                // Convert subkey to multistring.
                //
                newNode = (PININODE)malloc(sizeof(ININODE) + 
                                                  strlen(currentNode->Name) + 1 + 
                                                  multiLen + 1);
                if (newNode == NULL) {
                    status = STATUS_INI_INSUFFICIENT_RESOURCES;
                    goto fail;
                }
                newNode->Magic = ININODE_MAGIC;
                newNode->Type = INI_MULTI_SZ;
				newNode->LineNumber = keyLine;
                newNode->Name = (char*)(newNode + 1);
                newNode->ValueSize = multiLen + 1;
                newNode->Value = newNode->Name + strlen(currentNode->Name) + 1;
                newNode->Parent = currentNode->Parent;
				newNode->allocatedValue = false;
                strcpy(newNode->Name, currentNode->Name);
                memmove(newNode->Value, multi, multiLen);
                newNode->Value[multiLen] = 0;
				newNode->Next = NULL;
                for (node = parentNode->SubKey; node != NULL; node = node->Next) {
                    if (node->Next == currentNode) {
						node->Next = newNode;
                        break;
                    }
                }
				free(currentNode);

				currentNode = newNode;

                parseState = 0;
                continue;
            }

            IniParseQuotedString(line);
            memmove(multi + multiLen, line, strlen(line)+1);
            multiLen += strlen(multi + multiLen) + 1;
            line = "";

            break;

        }

        while (line != NULL && *line == '\0') {
            line = IniReadLine(iniFile);
            lineNumber++;
        }
    }

    if (parentNode != &iniFile->RootNode || parseState != 0) {
        //
        // Unexpected end of file found.
        //
        status = STATUS_INI_FILE_UNEXPECTED_EOF + (lineNumber & 0xffffff);
        goto fail;
    }

    free(iniFile->Buffer);
    iniFile->Buffer = NULL;

    *hIniFile = iniFile;

    return bSignature ? STATUS_INI_SIGNATURE_OK : (bChecksum ? STATUS_INI_CHECKSUM_OK : STATUS_INI_SUCCESS);

fail:
	if (pPubKeyID)
		free(pPubKeyID);
	if (pSignature)
		free(pSignature);
	if (hMD5 != NULL)
		CryptoDestroyHash(hMD5);
	if (hSHA1 != NULL)
		CryptoDestroyHash(hSHA1);
    if (iniFile) {
        IniReleaseFile(iniFile);
    } else {
        free(buffer);
    }

    return status;

} // IniParseFile



VOID
IniReleaseFile(
    HANDLE hIniFile
    )

/*++

Routine Description:

    This routine releases resources allocated by IniParseFile.

Arguments:

    hIniFile - handle returned by IniParseFile.

Return Value:

    None.

--*/

{
    PINIFILE iniFile = (PINIFILE)hIniFile;
    PININODE node;

    if (iniFile && iniFile->Magic == INIFILE_MAGIC) {
        free(iniFile->Buffer);

        node = &iniFile->RootNode;
        while (iniFile->RootNode.SubKey != NULL) {
            if (node->Type == INI_NONE && node->SubKey != NULL) {
                node = node->SubKey;
            } else {
                PININODE next;
                next = node->Next;
                node->Parent->SubKey = next;
                if (next == NULL) {
                    next = node->Parent;
                }
				if(node->allocatedValue)
					free(node->Value);
                free(node);
                node = next;
            }
        }

        iniFile->Magic = 0;
        free(iniFile);
    }

} // IniReleaseFile



INISTATUS
IniOpenKey(
    HANDLE  hKey,
    const char *pSubKey,
    HANDLE *phkResult
    )

/*++

Routine Description:

    This function opens the specified key.

Arguments:

    hKey      - Identifies a currently open key.
    pSubKey  - Pointer to a null-terminated string containing the name 
                of the subkey to open.
    phkResult - Pointer to a variable that receives the handle 
                of the opened key.

Return Value:

    STATUS_SUCCESS
    STATUS_INVALID_HANDLE
    STATUS_NOT_FOUND

--*/

{
    PININODE node, subKey;
    const char *keyName;

    if (hKey == NULL) {
        return STATUS_INI_INVALID_HANDLE;
    }
    if (((PINIFILE)hKey)->Magic == INIFILE_MAGIC) {
        node = ((PINIFILE)hKey)->RootNode.SubKey;
    } else if (((PININODE)hKey)->Magic == ININODE_MAGIC || ((PININODE)hKey)->Magic == INISECT_MAGIC) {
        node = ((PININODE)hKey)->SubKey;
    } else {
        return STATUS_INI_INVALID_HANDLE;
    }

    *phkResult = NULL;

    while (node != NULL) {
        for (subKey = node; subKey != NULL; subKey = subKey->Next) {
             if (strcmp(pSubKey, subKey->Name) == 0 && subKey->Type == INI_NONE) {
                 *phkResult = (HANDLE)subKey;
                 return STATUS_INI_SUCCESS;
             }
        }
        keyName = pSubKey;
        while (*pSubKey != '\0' && *pSubKey != '\\') {
            pSubKey++;
        }
        if (*pSubKey == '\0') {
            return STATUS_INI_KEY_NOT_FOUND;
        }
        for (subKey = node; subKey != NULL; subKey = subKey->Next) {
			if (strncmp(keyName, subKey->Name, pSubKey - keyName) == 0) {
				if (subKey->Name[pSubKey - keyName] == 0) {
					break;
				}
			}
        }
        if (subKey == NULL) {
            return STATUS_INI_KEY_NOT_FOUND;
        }
		pSubKey++;
    }

    return STATUS_INI_KEY_NOT_FOUND;

}



INISTATUS
IniEnumKey(
    HANDLE  hKey,
    DWORD   dwIndex,
    const char **ppName,
	DWORD  *pLineNumber
    )

/*++

Routine Description:

    This function enumerates subkeys of the specified key. The function 
    retrieves one subkey each time it is called. 

Arguments:

    hKey     - handle of key to enumerate returned by IniOpenKey or IniParseFile.
    dwIndex  - index of subkey to enumerate. Index should be zero for the first
               call and then incremented for subsequent calls. 
    ppName   - address of pointer to subkey name.

Return Value:

    STATUS_SUCCESS
    STATUS_INI_INVALID_HANDLE
    STATUS_NO_MORE_ENTRIES

--*/

{
    PININODE node, subKey;
    DWORD index;

    if (hKey == NULL) {
        return STATUS_INI_INVALID_HANDLE;
    }
    if (((PINIFILE)hKey)->Magic == INIFILE_MAGIC) {
        node = &(((PINIFILE)hKey)->RootNode);
    } else if (((PININODE)hKey)->Magic == ININODE_MAGIC || ((PININODE)hKey)->Magic == INISECT_MAGIC) {
        node = (PININODE)hKey;
    } else {
        return STATUS_INI_INVALID_HANDLE;
    }

    for (index = 0, subKey = node->SubKey; 
         subKey != NULL; 
         subKey = subKey->Next) 
    {
         if (subKey->Type == INI_NONE) {
             if (index == dwIndex) {
                 break;
             }
             index++;
         }
    }

    if (subKey == NULL) {
        return STATUS_INI_NO_MORE_ENTRIES;
    }

    *ppName = subKey->Name;
    if (pLineNumber != NULL) {
        *pLineNumber = subKey->LineNumber;
    }

    return STATUS_INI_SUCCESS;

} // IniEnumKey



INISTATUS 
IniEnumValue(
    HANDLE       hKey,
    DWORD        dwIndex,
    const char **ppName,
    DWORD        *pType,
    const char **ppData,
    DWORD       *pcbData,
	DWORD       *pLineNumber
    )

/*++

Routine Description:

    This function enumerates subkeys of the specified key. The function 
    retrieves one subkey each time it is called. 

Arguments:

    hKey     - handle of key to enumerate returned by IniOpenKey or IniParseFile.
    dwIndex  - index of subkey to enumerate. Index should be zero for the first
               call and then incremented for subsequent calls. 
    ppName   - pointer to a variable that receives a pointer to the subkey name.
    pType    - pointer to a variable that receives the type code for the value.
               The lpType can be NULL if the type code is not required. 
    ppData   - pointer to a variable that receives a pointer to the data 
	           of the value entry. Can be NULL if the data is not required. 
    pcbData  - pointer to a variable that receives the size in bytes 
	           of the data for the value entry. The pcbData can be NULL.

Return Value:

    STATUS_SUCCESS
    STATUS_INI_INVALID_HANDLE
    STATUS_NO_MORE_ENTRIES

--*/

{
    PININODE node, subKey;
    DWORD index;

    if (hKey == NULL) {
        return STATUS_INI_INVALID_HANDLE;
    }
    if (((PINIFILE)hKey)->Magic == INIFILE_MAGIC) {
        node = &(((PINIFILE)hKey)->RootNode);
    } else if (((PININODE)hKey)->Magic == ININODE_MAGIC || ((PININODE)hKey)->Magic == INISECT_MAGIC) {
        node = (PININODE)hKey;
    } else {
        return STATUS_INI_INVALID_HANDLE;
    }

    for (index = 0, subKey = node->SubKey; 
         subKey != NULL; 
         subKey = subKey->Next) 
    {
        if (subKey->Type != INI_NONE) {
            if (index == dwIndex) {
                break;
            }
            index++;
        }
    }

    if (subKey == NULL) {
        return STATUS_INI_NO_MORE_ENTRIES;
    }

    *ppName = subKey->Name;
    if (pType != NULL) {
        *pType = subKey->Type;
    }
    if (ppData != NULL) {
        *ppData = subKey->Value;
    }
    if (pcbData != NULL) {
        *pcbData = subKey->ValueSize;
    }
    if (pLineNumber != NULL) {
        *pLineNumber = subKey->LineNumber;
    }

    return STATUS_INI_SUCCESS;

} // IniEnumValue



INISTATUS
IniQueryInfoKey(
    HANDLE  hKey,
    DWORD  *pcSubKeys,
    DWORD  *pcValues,
	DWORD  *pLineNumber
    )

/*++

Routine Description:

    This function retrieves information about a specified registry key.

Arguments:

    hKey                - Identifies a currently open key.
    lpcSubKeys          - Address of buffer for number of  subkeys
    lpcValues           - Address of buffer for number of value entries.

Return Value:

    STATUS_SUCCESS
    STATUS_INI_INVALID_HANDLE

--*/

{
    PININODE node;

    if (hKey == NULL) {
        return STATUS_INI_INVALID_HANDLE;
    }
    if (((PINIFILE)hKey)->Magic == INIFILE_MAGIC) {
        node = ((PINIFILE)hKey)->RootNode.SubKey;
    } else if (((PININODE)hKey)->Magic == ININODE_MAGIC || ((PININODE)hKey)->Magic == INISECT_MAGIC) {
        node = ((PININODE)hKey)->SubKey;
    } else {
        return STATUS_INI_INVALID_HANDLE;
    }

    if (pcSubKeys != NULL) {
        *pcSubKeys = 0;
    }
    if (pcValues != NULL) {
        *pcValues = 0;
    }
    if (pLineNumber != NULL) {
        *pLineNumber = node->LineNumber;
    }

    while (node != NULL) {
        if (node->Type == INI_NONE) {
            if (pcSubKeys != NULL) {
                (*pcSubKeys)++;
            }
        } else {
            if (pcValues != NULL) {
                (*pcValues)++;
            }
        }
        node = node->Next;
    }

    return STATUS_INI_SUCCESS;

} // IniQueryInfoKey



INISTATUS 
IniQueryValue(
    HANDLE       hKey,
    const char  *pValueName,
    DWORD       *pType,
    const char **ppData,
    DWORD       *pcbData,
	DWORD		*pLineNumber
    )

/*++

Routine Description:

    This function enumerates subkeys of the specified key. The function 
    retrieves one subkey each time it is called. 

Arguments:

    hKey        - Identifies a currently open key.
    pValueName  - Pointer to a null-terminated string containing the name 
                  of the value to query.
    lpType      - pointer to a variable that receives the type code for 
                  the value. The lpType can be NULL if the type is not required. 
    ppData      - pointer to a buffer that receives the data for the value 
                  entry. Can be NULL if the data is not required. 
    pcbData     - pointer to a variable that receives the size in bytes 
	              of the data for the value entry. The pcbData can be NULL.

Return Value:

    STATUS_INI_SUCCESS
    STATUS_INI_INVALID_HANDLE
    STATUS_INI_BUFFER_OVERFLOW
    STATUS_INI_KEY_NOT_FOUND

--*/

{
    PININODE node, subKey;

    if (hKey == NULL) {
        return STATUS_INI_INVALID_HANDLE;
    }
    if (((PINIFILE)hKey)->Magic == INIFILE_MAGIC) {
        node = &((PINIFILE)hKey)->RootNode;
    } else if (((PININODE)hKey)->Magic == ININODE_MAGIC || ((PININODE)hKey)->Magic == INISECT_MAGIC) {
        node = (PININODE)hKey;
    } else {
        return STATUS_INI_INVALID_HANDLE;
    }

    for (subKey = node->SubKey; subKey != NULL; subKey = subKey->Next) {
        if (subKey->Type != INI_NONE && strcmp(subKey->Name, pValueName) == 0)
            break;
    }

    if (subKey == NULL) {
        return STATUS_INI_KEY_NOT_FOUND;
    }

    if (pType != NULL) {
        *pType = subKey->Type;
    }
    if (ppData != NULL) {
        *ppData = subKey->Value;
    }
    if (pcbData != NULL) {
        *pcbData = subKey->ValueSize;
	}
    if (pLineNumber != NULL) {
        *pLineNumber = subKey->LineNumber;
    }

    return STATUS_INI_SUCCESS;

} // IniQueryValue


const char *
IniErrorMessage(
	INISTATUS status
	)

/*++

Routine Description:

    This function return error message coreresponding to error code

Arguments:

    status      - Error code.

Return Value:

    Error mesage.

--*/

{
	if (INI_SYNTAX_ERROR(status)) {
		DWORD line = INI_ERROR_LINE_NUMBER(status);
		status -= line;
		switch (status) {
		case STATUS_INI_FILE_NOT_TOP_LEVEL:
			return "section header not at the top level";
		case STATUS_INI_FILE_KEY_EXISTS:
			return "duplicate key found";
		case STATUS_INI_FILE_EXTRA_CBRACE:
			return "extra closing brace";
		case STATUS_INI_FILE_NEWLINE:
			return "missing closing quotation mark";
		case STATUS_INI_FILE_UNEXPECTED_EOF:
			return "unexpected end of file";
		case STATUS_INI_FILE_SYNTAX_ERROR:
		default:
			return "syntax error";
		}
	} else {
		switch (status) {
		case STATUS_INI_FILE_READ_FAILED:
			return "file read failed";
		case STATUS_INI_INSUFFICIENT_RESOURCES:
			return "not enough memory";
		case STATUS_INI_KEY_NOT_FOUND:
			return "key not found";
		case STATUS_INI_FILE_IS_TOO_BIG:
			return "file is too big";
		case STATUS_INI_FILE_NOT_FOUND:
			return "file not found";
		default:
			return "";
		}
	}
}



INISTATUS 
IniCreateFile(
	HANDLE* phIniFile
	)

/*++

Routine Description:

    This routine creates a new ini file.

Arguments:

	phIniFile - Pointer to a variable that receives a handle to the created init file.

Return Value:

    STATUS_INI_SUCCESS
    STATUS_INI_INSUFFICIENT_RESOURCES

--*/

{
	INIFILE* iniFile = (PINIFILE)malloc(sizeof(INIFILE));

	if (iniFile == NULL) {
        return STATUS_INI_INSUFFICIENT_RESOURCES;
    }

    iniFile->Magic = INIFILE_MAGIC;
    iniFile->Buffer = NULL;
    iniFile->BufferSize = 0;
    iniFile->BufferPos = 0;
    iniFile->RootNode.Magic = ININODE_MAGIC;
    iniFile->RootNode.Type = INI_NONE;
    iniFile->RootNode.Next = NULL;
    iniFile->RootNode.Parent = NULL;
    iniFile->RootNode.Name = NULL;
    iniFile->RootNode.ValueSize = 0;
    iniFile->RootNode.Value = NULL;
	iniFile->RootNode.allocatedValue = false;
	iniFile->RootNode.SubKey = NULL;

	if (phIniFile)
		*phIniFile = (HANDLE)iniFile;

	return STATUS_INI_SUCCESS;
}



INISTATUS
IniCreateKey(
	HANDLE	   hKey,
	const char *pKeyName,
	HANDLE	  *phkResult
	)

/*++

Routine Description:

    This routine creates a new key. If the key already exists routine opens it.

Arguments:

    hKey - Identifies a currently open key
	pKeyName - Name of the key to create.
	phkResult - Pointer to a variable that receives a handle to the created key.

Return Value:

    STATUS_INI_SUCCESS
    STATUS_INI_INVALID_HANDLE
    STATUS_INI_INSUFFICIENT_RESOURCES

--*/

{
	ININODE *parentNode;
	ININODE *newNode;
	ININODE *node;
	
    if (hKey == NULL)
        return STATUS_INI_INVALID_HANDLE;
    if (((PINIFILE)hKey)->Magic == INIFILE_MAGIC) {
        parentNode = &((PINIFILE)hKey)->RootNode;
    } else if (((PININODE)hKey)->Magic == ININODE_MAGIC || ((PININODE)hKey)->Magic == INISECT_MAGIC) {
        parentNode = (PININODE)hKey;
	} else {
        return STATUS_INI_INVALID_HANDLE;
	}

    for (node = parentNode->SubKey; node != NULL; node = node->Next) {
        if (strcmp(node->Name, pKeyName) == 0) {
			if (phkResult)
				*phkResult = (HANDLE)node;
            return STATUS_INI_SUCCESS;
        }
    }

    //
    // Add new key
    //
    newNode = (PININODE)malloc(sizeof(ININODE) + strlen(pKeyName) + 1);
    if (newNode == NULL) {
        return STATUS_INI_INSUFFICIENT_RESOURCES;
    }
    newNode->Magic = ININODE_MAGIC;
    newNode->Type = INI_NONE;
	newNode->LineNumber = 0;
    newNode->Name = (char*)(newNode + 1);
    newNode->Value = NULL;
    newNode->Parent = parentNode;
	newNode->allocatedValue = false;
    strcpy(newNode->Name, pKeyName);
    for (node = parentNode->SubKey; node != NULL; node = node->Next) {
        if (node->Next == NULL) {
            break;
        }
    }
    newNode->Next = NULL;
    if (node == NULL) {
        parentNode->SubKey = newNode;
    } else {
        node->Next = newNode;
    }

	if (phkResult)
		*phkResult = (HANDLE)newNode;

	return STATUS_INI_SUCCESS;
}



INISTATUS
IniSetValue(
	HANDLE	    hKey,
	const char *pValueName,
	DWORD       dwType,
	const char *lpData
	)

/*++

Routine Description:

    This routine sets the data and type of a specified value under a key.

Arguments:

    hKey - Identifies a currently open key
	pValueName - Name of the value to set.
	dwType - Specifies a code indicating the type of value
	lpData - Data to be stored with the specified value name.

Return Value:

    STATUS_INI_SUCCESS
    STATUS_INI_INVALID_HANDLE
    STATUS_INI_INVALID_VALUE_TYPE
    STATUS_INI_INSUFFICIENT_RESOURCES

--*/

{
	ININODE *parentNode;
	ININODE *newNode;
	ININODE *node;
	
    if (hKey == NULL)
        return STATUS_INI_INVALID_HANDLE;
    if (((PINIFILE)hKey)->Magic == INIFILE_MAGIC) {
        parentNode = &((PINIFILE)hKey)->RootNode;
    } else if (((PININODE)hKey)->Magic == ININODE_MAGIC || ((PININODE)hKey)->Magic == INISECT_MAGIC) {
        parentNode = (PININODE)hKey;
	} else {
        return STATUS_INI_INVALID_HANDLE;
	}

	if (dwType != INI_SZ && dwType != INI_MULTI_SZ)
		return STATUS_INI_INVALID_VALUE_TYPE;

    for (node = parentNode->SubKey; node != NULL; node = node->Next) {
        if (node->Type != INI_NONE && strcmp(node->Name, pValueName) == 0)
            break;
    }

	DWORD cbData = 0;
	if (dwType == INI_SZ) {
		cbData = strlen(lpData) + 1;
	} else if (dwType == INI_MULTI_SZ) {
		const char *p = lpData;
		while (*p || *(p+1)) {
			int len = strlen(p) + 1;
			cbData += len;
			p += len;
			if (*p == 0)
				break;
		}
		cbData += 1;
	}

	if (node != NULL) {
		//
		// Update value
		//
		char *newValue = (char*)malloc(cbData);
		if (newValue == NULL) {
			return STATUS_INI_INSUFFICIENT_RESOURCES;
		}
		memmove(newValue, lpData, cbData);
		if (node->allocatedValue)
			free(node->Value);
		node->Value = newValue;
		node->allocatedValue = true;
	} else {
		//
		// Add new value
		//
		newNode = (PININODE)malloc(sizeof(ININODE) + strlen(pValueName) + 1);
		if (newNode == NULL) {
			return STATUS_INI_INSUFFICIENT_RESOURCES;
		}
		newNode->Value = (char*)malloc(cbData);
		if (newNode->Value == NULL) {
			free(newNode);
			return STATUS_INI_INSUFFICIENT_RESOURCES;
		}
		memmove(newNode->Value, lpData, cbData);
		newNode->allocatedValue = true;
		newNode->Magic = ININODE_MAGIC;
		newNode->Type = dwType;
		newNode->LineNumber = 0;
		newNode->Name = (char*)(newNode + 1);
		newNode->Parent = parentNode;
		strcpy(newNode->Name, pValueName);
		for (node = parentNode->SubKey; node != NULL; node = node->Next) {
			if (node->Next == NULL) {
				break;
			}
		}
		newNode->Next = NULL;
		if (node == NULL) {
			parentNode->SubKey = newNode;
		} else {
			node->Next = newNode;
		}
	}

	return STATUS_INI_SUCCESS;

} // IniSetValue



INISTATUS
IniDeleteKey(
    HANDLE hKey,
	const char *pSubKey
    )

/*++

Routine Description:

    This routine deletes a key and all all its subkey.

Arguments:

    hKey - Identifies a currently open key.
	pSubKey - Name of a subkey to delete.

Return Value:

    STATUS_INI_SUCCESS
    STATUS_INI_INVALID_HANDLE

--*/

{
	ININODE *parentNode;
	ININODE *keyNode;
	ININODE *node;
	
    if (hKey == NULL)
        return STATUS_INI_INVALID_HANDLE;
    if (((PINIFILE)hKey)->Magic == INIFILE_MAGIC) {
        parentNode = &((PINIFILE)hKey)->RootNode;
    } else if (((PININODE)hKey)->Magic == ININODE_MAGIC || ((PININODE)hKey)->Magic == INISECT_MAGIC) {
        parentNode = (PININODE)hKey;
	} else {
        return STATUS_INI_INVALID_HANDLE;
	}

    node = parentNode->SubKey;
	if (node == NULL) {
		return STATUS_INI_KEY_NOT_FOUND;
	} if (node->Type == INI_NONE && strcmp(node->Name, pSubKey) == 0) {
		parentNode->SubKey = node->Next;
	} else {
		while (node != NULL) {
			if (node->Next 
			 && node->Next->Type == INI_NONE 
			 && strcmp(node->Next->Name, pSubKey) == 0) 
			{
				node->Next = node->Next->Next;
				node = node->Next;
				break;
			}
			node = node->Next;
		}
		if (node == NULL)
			return STATUS_INI_KEY_NOT_FOUND;
	}

	keyNode = node;
    node = keyNode->SubKey;
    while (keyNode->SubKey != NULL) {
        if (node->Type == INI_NONE && node->SubKey != NULL) {
            node = node->SubKey;
        } else {
            PININODE next;
            next = node->Next;
            node->Parent->SubKey = next;
            if (next == NULL) {
                next = node->Parent;
            }
			if(node->allocatedValue)
				free(node->Value);
            free(node);
            node = next;
        }
    }
	if(keyNode->allocatedValue)
		free(keyNode->Value);
	free(keyNode);

	return STATUS_INI_SUCCESS;

} // IniDeleteKey



INISTATUS
IniDeleteValue(
    HANDLE hKey,
	const char *pValueName
    )

/*++

Routine Description:

    This routine deletes a named value from the specified key.

Arguments:

    hKey - Identifies a currently open key.
	pValueName - Name of a value to delete.

Return Value:

    STATUS_INI_SUCCESS
    STATUS_INI_INVALID_HANDLE

--*/

{
	ININODE *parentNode;
	ININODE *node;
	
    if (hKey == NULL)
        return STATUS_INI_INVALID_HANDLE;
    if (((PINIFILE)hKey)->Magic == INIFILE_MAGIC) {
        parentNode = &((PINIFILE)hKey)->RootNode;
    } else if (((PININODE)hKey)->Magic == ININODE_MAGIC || ((PININODE)hKey)->Magic == INISECT_MAGIC) {
        parentNode = (PININODE)hKey;
	} else {
        return STATUS_INI_INVALID_HANDLE;
	}

    node = parentNode->SubKey;
	if (node == NULL) {
		return STATUS_INI_KEY_NOT_FOUND;
	} if (node->Type != INI_NONE && strcmp(node->Name, pValueName) == 0) {
		parentNode->SubKey = node->Next;
	} else {
		while (node != NULL) {
			if (node->Next 
			 && node->Next->Type != INI_NONE 
			 && strcmp(node->Next->Name, pValueName) == 0) 
			{
				node->Next = node->Next->Next;
				node = node->Next;
				break;
			}
			node = node->Next;
		}
		if (node == NULL)
			return STATUS_INI_KEY_NOT_FOUND;
	}
	if(node->allocatedValue)
		free(node->Value);
	free(node);

	return STATUS_INI_SUCCESS;

} // IniDeleteValue


static bool IniWriteBuf(BUF_STATE *pState, const char *pBuffer, int nBufferSize)
{
	if (pState->stringBuffer) {
		pState->stringBuffer->append(pBuffer, nBufferSize);
	} else if (write(pState->fileHandle, pBuffer, nBufferSize) != nBufferSize) {
		return false;
	}
	if (pState->hSHA1)
		CryptoHashData(pState->hSHA1, (BYTE*)pBuffer, nBufferSize);
	if (pState->hMD5)
		CryptoHashData(pState->hMD5, (BYTE*)pBuffer, nBufferSize);
	return true;
}



INISTATUS
IniWrite(
    HANDLE hIniFile,
    BUF_STATE *bufState,
	DWORD dwChecksumType,
	INI_MAKE_SIGNATURE_CALLBACK pSignCallback,
	void *pSignCallbackCtx
    )

/*++

Routine Description:

	This routine writes ini file.

Arguments:

    hIniFile - handle returned by IniParseFile of IniCreateFile.
	pFileName - name of of the file to write to.

Return Value:

    STATUS_INI_SUCCESS
    STATUS_INI_INVALID_HANDLE
	STATUS_INI_FILE_WRITE_FAILED

--*/

{
	INISTATUS status;
	INIFILE* pIniFile;
	BYTE hash[16+20];
	char hashhex[2*(16+20)+1];
	const char *pszHashType;
	DWORD nHashLen;
	DWORD dwLen;

    if (hIniFile == NULL || ((PINIFILE)hIniFile)->Magic != INIFILE_MAGIC)
        return STATUS_INI_INVALID_HANDLE;
	pIniFile = (INIFILE*)hIniFile;

	if (dwChecksumType == INI_CHECKSUM_MD5) {
		if (!CryptoCreateHash(CRYPTO_ALG_MD5, &bufState->hMD5))
			return STATUS_INI_INSUFFICIENT_RESOURCES;
		nHashLen = 16;
		pszHashType = CHECKSUM_HASH_MD5;
	} else if (dwChecksumType == INI_CHECKSUM_SHA1) {
		if (!CryptoCreateHash(CRYPTO_ALG_SHA1, &bufState->hSHA1))
			return STATUS_INI_INSUFFICIENT_RESOURCES;
		nHashLen = 20;
		pszHashType = CHECKSUM_HASH_SHA1;
	} else if (dwChecksumType == INI_CHECKSUM_MD5_SHA1) {
		if (!CryptoCreateHash(CRYPTO_ALG_MD5, &bufState->hMD5))
			return STATUS_INI_INSUFFICIENT_RESOURCES;
		if (!CryptoCreateHash(CRYPTO_ALG_SHA1, &bufState->hSHA1))
			return STATUS_INI_INSUFFICIENT_RESOURCES;
		nHashLen = 16+20;
		pszHashType = CHECKSUM_HASH_MD5_SHA1;
	} else if (dwChecksumType == INI_CHECKSUM_NONE) {
		nHashLen = 0;
		pszHashType = NULL;
	} else {
		return STATUS_INI_UNSUPPORTED_CHECKSUM;
	}

	status = IniWriteNode(bufState, pIniFile->RootNode.SubKey, 0);
	if (INI_FAILED(status))
		return status;

	if (dwChecksumType != INI_CHECKSUM_NONE) {
		DWORD dwHashLen;
		if (dwChecksumType == INI_CHECKSUM_MD5 || dwChecksumType == INI_CHECKSUM_MD5_SHA1) {
			dwHashLen = 16;
			if (!CryptoGetHashParam(bufState->hMD5, CRYPTO_HP_HASHVALUE, (BYTE*)hash, &dwHashLen)) {
				CryptoDestroyHash(bufState->hMD5);
				return STATUS_INI_INSUFFICIENT_RESOURCES;
			}
			CryptoDestroyHash(bufState->hMD5);
			bufState->hMD5 = NULL;
		}
		if (dwChecksumType == INI_CHECKSUM_SHA1 || dwChecksumType == INI_CHECKSUM_MD5_SHA1) {
			dwHashLen = 20;
			if (!CryptoGetHashParam(bufState->hSHA1, CRYPTO_HP_HASHVALUE, 
				(BYTE*)hash + (dwChecksumType == INI_CHECKSUM_MD5_SHA1 ? 16 : 0), &dwHashLen)) {
				CryptoDestroyHash(bufState->hSHA1);
				return STATUS_INI_INSUFFICIENT_RESOURCES;
			}
			CryptoDestroyHash(bufState->hSHA1);
			bufState->hSHA1 = NULL;
		}
		dwLen = sizeof(hashhex);
		if (!CryptoConvertToHex((BYTE*)hash, nHashLen, hashhex, &dwLen)) {
			return STATUS_INI_INSUFFICIENT_RESOURCES;
		}

		if (!IniWriteBuf(bufState, CHECKSUM_BEGIN, strlen(CHECKSUM_BEGIN))) {
			return STATUS_INI_FILE_WRITE_FAILED;
		}
		if (!IniWriteBuf(bufState, pszHashType, strlen(pszHashType))) {
			return STATUS_INI_FILE_WRITE_FAILED;
		}
		if (!IniWriteBuf(bufState, hashhex, strlen(hashhex))) {
			return STATUS_INI_FILE_WRITE_FAILED;
		}
		if (!IniWriteBuf(bufState, CHECKSUM_END, strlen(CHECKSUM_END))) {
			return STATUS_INI_FILE_WRITE_FAILED;
		}
	}

	if (pSignCallback) {
		DWORD dwSigType;
		BYTE *pPubKeyID;
		DWORD cbPubKeyID;
		BYTE *pSignature;
		DWORD cbSignature;
		const char *pszSigType;
		char *pszSigText;
		DWORD dwSigTextLen;
		
		status = pSignCallback(pSignCallbackCtx, dwChecksumType, hash, nHashLen,
			&dwSigType, &pPubKeyID, &cbPubKeyID, &pSignature, &cbSignature);
		if (INI_FAILED(status))
			return status;

		if (dwSigType == INI_SIGNATURE_NONE) {
			free(pPubKeyID);
			free(pSignature);
			return status;
		} else if (dwSigType == INI_SIGNATURE_RSA) {
			pszSigType = SIGNATURE_RSA;
		} else if (dwSigType == INI_SIGNATURE_DSA) {
			pszSigType = SIGNATURE_DSA;
		} else {
			free(pPubKeyID);
			free(pSignature);
			return STATUS_INI_UNSUPPORTED_SIGNATURE;
		}

		dwSigTextLen = 1;
		if (pPubKeyID) {
			if (!CryptoConvertToBase64(pPubKeyID, cbPubKeyID, NULL, &dwLen)) {
				free(pPubKeyID);
				free(pSignature);
				return STATUS_INI_INSUFFICIENT_RESOURCES;
			}
			dwSigTextLen += dwLen + 2;
		}
		if (!CryptoConvertToBase64(pSignature, cbSignature, NULL, &dwLen)) {
			free(pPubKeyID);
			free(pSignature);
			return STATUS_INI_INSUFFICIENT_RESOURCES;
		}
		dwSigTextLen += dwLen;

		pszSigText = (char*)malloc(dwSigTextLen);
		if (pszSigText == NULL) {
			free(pPubKeyID);
			free(pSignature);
			return STATUS_INI_INSUFFICIENT_RESOURCES;
		}

		if (pPubKeyID) {
			*pszSigText = '[';
			dwLen = dwSigTextLen - 2;
			if (!CryptoConvertToBase64(pPubKeyID, cbPubKeyID, pszSigText+1, &dwLen)) {
				free(pPubKeyID);
				free(pSignature);
				return STATUS_INI_INSUFFICIENT_RESOURCES;
			}
			free(pPubKeyID);
			strcat(pszSigText, "]");
		} else {
			*pszSigText = '\0';
		}
		dwLen = dwSigTextLen - strlen(pszSigText);
		if (!CryptoConvertToBase64(pSignature, cbSignature, pszSigText+strlen(pszSigText), &dwLen)) {
			free(pSignature);
			return STATUS_INI_INSUFFICIENT_RESOURCES;
		}
		free(pSignature);

		if (!IniWriteBuf(bufState, SIGNATURE_BEGIN, strlen(SIGNATURE_BEGIN))) {
			free(pszSigText);
			return STATUS_INI_FILE_WRITE_FAILED;
		}
		if (!IniWriteBuf(bufState, pszSigType, strlen(pszSigType))) {
			return STATUS_INI_FILE_WRITE_FAILED;
		}
		if (!IniWriteBuf(bufState, pszSigText, strlen(pszSigText))) {
			return STATUS_INI_FILE_WRITE_FAILED;
		}
		if (!IniWriteBuf(bufState, SIGNATURE_END, strlen(SIGNATURE_END))) {
			return STATUS_INI_FILE_WRITE_FAILED;
		}
		
	}

	return status;

} // IniWrite



INISTATUS
IniWriteFile(
    HANDLE hIniFile,
    const char *pFileName,
	DWORD dwChecksumType,
	INI_MAKE_SIGNATURE_CALLBACK pSignCallback,
	void *pSignCallBackCtx
    )

/*++

Routine Description:

	This routine writes ini file.

Arguments:

    hIniFile - handle returned by IniParseFile of IniCreateFile.
	pFileName - name of of the file to write to.

Return Value:

    STATUS_INI_SUCCESS
    STATUS_INI_INVALID_HANDLE
	STATUS_INI_FILE_WRITE_FAILED

--*/

{
	INISTATUS status;
	INIFILE* pIniFile;
	BUF_STATE bufState;

    if (hIniFile == NULL || ((PINIFILE)hIniFile)->Magic != INIFILE_MAGIC)
        return STATUS_INI_INVALID_HANDLE;
	pIniFile = (INIFILE*)hIniFile;

#ifdef _WIN32
	int file = open(pFileName, O_RDWR | O_BINARY | O_TRUNC);
	if (file == -1)
		file = open(pFileName, O_RDWR | O_BINARY | O_CREAT, S_IREAD | S_IWRITE);	
#else
	int file = open(pFileName, O_RDWR | O_TRUNC);
	if (file == -1)
		file = open(pFileName, O_RDWR | O_CREAT, S_IREAD | S_IWRITE);	
#endif
	if (file == -1)
		return STATUS_INI_INVALID_HANDLE;

	bufState.stringBuffer = NULL;
	bufState.fileHandle = file;
	bufState.hMD5  = NULL;
	bufState.hSHA1 = NULL;

	status = IniWrite(hIniFile, &bufState, dwChecksumType, pSignCallback, pSignCallBackCtx);
	close(bufState.fileHandle);
	if (status != STATUS_INI_SUCCESS)
		unlink(pFileName);

	if (bufState.hMD5)
		CryptoDestroyHash(bufState.hMD5);
	if (bufState.hSHA1)
		CryptoDestroyHash(bufState.hSHA1);

	return status;
}



INISTATUS
IniWriteBuffer(
    HANDLE hIniFile,
    string &strBuffer,
	DWORD dwChecksumType,
	INI_MAKE_SIGNATURE_CALLBACK pSignCallback,
	void *pSignCallBackCtx
    )

/*++

Routine Description:

	This routine writes ini file.

Arguments:

    hIniFile - handle returned by IniParseFile of IniCreateFile.
	pFileName - name of of the file to write to.

Return Value:

    STATUS_INI_SUCCESS
    STATUS_INI_INVALID_HANDLE
	STATUS_INI_FILE_WRITE_FAILED

--*/

{
	INISTATUS status;
	INIFILE* pIniFile;
	BUF_STATE bufState;

    if (hIniFile == NULL || ((PINIFILE)hIniFile)->Magic != INIFILE_MAGIC)
        return STATUS_INI_INVALID_HANDLE;
	pIniFile = (INIFILE*)hIniFile;

	bufState.stringBuffer = &strBuffer;
	bufState.fileHandle = -1;
	bufState.hMD5  = NULL;
	bufState.hSHA1 = NULL;

	status = IniWrite(hIniFile, &bufState, dwChecksumType, pSignCallback, pSignCallBackCtx);
	if (status != STATUS_INI_SUCCESS)
		strBuffer = "";

	if (bufState.hMD5)
		CryptoDestroyHash(bufState.hMD5);
	if (bufState.hSHA1)
		CryptoDestroyHash(bufState.hSHA1);

	return status;
}


INISTATUS
IniWriteNode(
	BUF_STATE *pBufState,
	ININODE *pNode,
	int indentLevel)
{
	INISTATUS status;

#define INI_WRITE(buf, size) \
	{ \
        if (!IniWriteBuf(pBufState, buf, size)) \
		    return STATUS_INI_FILE_WRITE_FAILED; \
	}

#define INI_WRITE_INDENT(level) \
	{ \
        for (int i = 0; i < level; i++) { \
			if (!IniWriteBuf(pBufState, "    ", 4)) \
    		    return STATUS_INI_FILE_WRITE_FAILED; \
		} \
	}

#define INI_WRITE_STRING(str) \
	{ \
		if (!IniWriteBuf(pBufState, "\"", 1)) \
    		return STATUS_INI_FILE_WRITE_FAILED; \
		char *p = str; \
		while (*p) { \
			switch (*p) { \
			case '"': \
				if (!IniWriteBuf(pBufState, "\\\"", 2)) \
					return STATUS_INI_FILE_WRITE_FAILED; \
				break; \
			case '\a': \
				if (IniWriteBuf(pBufState, "\\a", 2)) \
					return STATUS_INI_FILE_WRITE_FAILED; \
				break; \
			case '\b': \
				if (!IniWriteBuf(pBufState, "\\b", 2)) \
					return STATUS_INI_FILE_WRITE_FAILED; \
				break; \
			case '\f': \
				if (!IniWriteBuf(pBufState, "\\f", 2)) \
					return STATUS_INI_FILE_WRITE_FAILED; \
				break; \
			case '\n': \
				if (!IniWriteBuf(pBufState, "\\n", 2)) \
					return STATUS_INI_FILE_WRITE_FAILED; \
				break; \
			case '\r': \
				if (!IniWriteBuf(pBufState, "\\r", 2)) \
					return STATUS_INI_FILE_WRITE_FAILED; \
				break; \
			case '\t': \
				if (!IniWriteBuf(pBufState, "\\t", 2)) \
					return STATUS_INI_FILE_WRITE_FAILED; \
				break; \
			case '\v': \
				if (!IniWriteBuf(pBufState, "\\v", 2)) \
					return STATUS_INI_FILE_WRITE_FAILED; \
				break; \
			case '\\': \
				if (!IniWriteBuf(pBufState, "\\\\", 2)) \
					return STATUS_INI_FILE_WRITE_FAILED; \
				break; \
			default: \
				if (!IniWriteBuf(pBufState, p, 1)) \
					return STATUS_INI_FILE_WRITE_FAILED; \
				break; \
			} \
			p++; \
		} \
		if (!IniWriteBuf(pBufState, "\"", 1)) \
    		return STATUS_INI_FILE_WRITE_FAILED; \
	}

	for (; pNode != NULL; pNode = pNode->Next) {

		if (pNode->Type == INI_NONE) {

			if (indentLevel == 0) {
				INI_WRITE("[", 1)
				INI_WRITE(pNode->Name, strlen(pNode->Name))
				INI_WRITE("]\n", 2)
			} else {
				INI_WRITE_INDENT(indentLevel)
				INI_WRITE(pNode->Name, strlen(pNode->Name))
				INI_WRITE(" = {\n", 5)
			}
		
			status = IniWriteNode(pBufState, pNode->SubKey, indentLevel + 1);
			if (status != STATUS_INI_SUCCESS)
				return status;

			if (indentLevel != 0) {
				INI_WRITE_INDENT(indentLevel)
				INI_WRITE("}\n", 2);
			}

		} else if (pNode->Type == INI_SZ) {

			INI_WRITE_INDENT(indentLevel)
			INI_WRITE(pNode->Name, strlen(pNode->Name))
			INI_WRITE(" = ", 3)
			INI_WRITE_STRING(pNode->Value)
			INI_WRITE("\n", 1)

		} else if (pNode->Type == INI_MULTI_SZ) {

			INI_WRITE_INDENT(indentLevel)
			INI_WRITE(pNode->Name, strlen(pNode->Name))
			INI_WRITE(" = {\n", 5);

			char *pValue = pNode->Value;
			while (pValue) {
				INI_WRITE_INDENT(indentLevel + 1)
				INI_WRITE_STRING(pValue)
				INI_WRITE("\n", 1)
				pValue += strlen(pValue) + 1;
				if (*pValue == 0)
					break;
			}

			INI_WRITE_INDENT(indentLevel)
			INI_WRITE("}\n", 2)

		}

	}
	
	return STATUS_INI_SUCCESS;

#undef INI_WRITE
#undef INI_WRITE_INDENT
#undef INI_WRITE_STRING
}	

