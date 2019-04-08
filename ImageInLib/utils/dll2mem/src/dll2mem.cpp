#include <SDKDDKVer.h>
#include <Windows.h>
#include <stdint.h>
#include <stdio.h>
#include <io.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <share.h>
#include <time.h>
#include <string>
#include <vector>
#include <list>
#include <set>
#include <map>
#include <random>
#include "RLEEncoder.h"
#include "DiBitsEncoder.h"
#include "lzma/LzmaEnc.h"
#include "lzma/LzmaDec.h"


bool opt_SingleSection = false;
bool opt_LZMA = false;
bool opt_DiBits = true;
bool opt_64bit = false;
bool opt_Asm = false;

#define PAGE_SIZE 4096

struct SECTION_TABLE_ENTRY
{
	DWORD dwSectionBase;
	DWORD dwSectionSize;
	DWORD dwProtection;
	DWORD dwDataSize;
	std::string Data;
};

struct RELOCATION_BLOCK
{
	DWORD dwRelocationBase;
	DWORD dwRelocationEntries;
	BYTE  Relocations[64];
};

struct IMPORT_TABLE_ENTRY
{
	std::string ModuleName;
	std::string ImportName;
	DWORD dwImportOrdinal;
	std::string SymbolName;
};

struct IMPORT_BLOCK
{
	DWORD dwImportReference;
	std::list<IMPORT_TABLE_ENTRY> Imports;
};

struct EXPORT_TABLE_ENTRY
{
	DWORD dwAddressOfFunction;
	std::list<std::string> FunctionNames;
	DWORD dwExportOrdinal;
};

struct DLL
{
	bool is64bit;
	PBYTE pFileData;
	DWORD dwFileSize;
	PIMAGE_DOS_HEADER pDosHeader;
	PIMAGE_NT_HEADERS32 pNtHeaders32;
	PIMAGE_NT_HEADERS64 pNtHeaders64;
	DWORD dwAllocationSize;
	PBYTE pImageBase;
	std::vector<SECTION_TABLE_ENTRY> Sections;
	std::string Relocations;
	std::list<IMPORT_BLOCK> Imports;
	std::list<IMPORT_BLOCK> DllImports;
	std::vector<EXPORT_TABLE_ENTRY> Exports;

	DLL()
		: is64bit(false)
		, pFileData(nullptr)
		, dwFileSize(0)
		, pDosHeader(nullptr)
		, pNtHeaders32(nullptr)
		, pNtHeaders64(nullptr)
		, dwAllocationSize(0)
		, pImageBase(nullptr)
	{ }
};


struct COFF_SECTION
{
	std::string name;
	DWORD virtualSize;
	DWORD virtualAddress;
	DWORD rawDataSize;
	PBYTE rawData;
	DWORD relocationsCount;
	PIMAGE_RELOCATION relocations;
	DWORD characteristics;

	COFF_SECTION()
		: virtualSize(0)
		, virtualAddress(0)
		, rawDataSize(0)
		, rawData(nullptr)
		, relocationsCount(0)
		, relocations(nullptr)
		, characteristics(0)
	{ }
};

struct COFF_SYMBOL
{
	std::string name;
    DWORD   value;
    SHORT   sectionNumber;
    WORD    type;
    BYTE    storageClass;
	DWORD   auxDataSize;
	PIMAGE_AUX_SYMBOL auxData;

	COFF_SYMBOL() 
		: value(0)
		, sectionNumber(0)
		, type(0)
		, storageClass(0)
		, auxDataSize(0)
		, auxData(nullptr)
	{ }
};

struct COFF
{
	PBYTE pFileData;
	DWORD dwFileSize;
	bool is64bit;
	std::vector<COFF_SECTION> Sections;
	std::vector<COFF_SYMBOL> Symbols;

	COFF()
		: pFileData(nullptr)
		, dwFileSize(0)
	{ }
};

struct IMPORT_TRANSLATION_TABLE
{
	std::map<std::string, std::string> dllImports;
	std::map<std::string, std::string> imports;
};


static void *LzmaAlloc(void *p, size_t size)
{
	p = p; 
	if (size == 0)
		return 0;
	return malloc(size);
}

static void LzmaFree(void *p, void *address)
{
	p = p; 
	free(address);
}

static ISzAlloc g_Alloc = { LzmaAlloc, LzmaFree };

unsigned int GetLSFRBits(int nBits, unsigned int& lfsr)
{
	unsigned int o = 0;
	unsigned int l = lfsr;
	while (nBits--)
	{
		l = (l >> 1) ^ (unsigned int)(0 - (l & 1u) & 0xd0000001u);
		o = (o << 1) + (l & 1);
	}
	lfsr = l;
	return o;
}

unsigned int UpdateLCG(unsigned int& lcg, unsigned int a, unsigned c)
{
	lcg = a * lcg + c;
	return lcg;
}

#define LCG_134775813_1(lcg) UpdateLCG(lcg, 134775813, 1)
#define LCG_1103515245_12345(lcg) UpdateLCG(lcg, 1103515245, 12345)
#define LCG_22695477_1(lcg) UpdateLCG(lcg, 22695477, 1)
#define LCG_214013_2531011(lcg) UpdateLCG(lcg, 214013, 2531011)
#define LCG_69069_1(lcg) UpdateLCG(lcg, 69069, 1)
#define LCG_1664525_1013904223(lcg) UpdateLCG(lcg, 1664525, 1013904223)


DLL* LoadDLLToMemory(const char* pszDllFilePath)
{
	int fd;
	size_t size;
	PBYTE data = NULL;
	
	if (_sopen_s(&fd, pszDllFilePath, _O_BINARY | _O_RDONLY, _SH_DENYNO, 0))
	{
		fprintf(stderr, "error: could not open DLL file \"%s\"\n", pszDllFilePath);
		return NULL;
	}

	size = _filelength(fd);
	data = new BYTE[size];
	if (data == NULL)
	{
		_close(fd);
		fprintf(stderr, "error: could not allocate buffer for file \"%s\"\n", pszDllFilePath);
		return NULL;
	}
	if (_read(fd, data, size) != size)
	{
		delete data;
		_close(fd);
		fprintf(stderr, "error: could not read file \"%s\"\n", pszDllFilePath);
		return NULL;
	}
	_close(fd);

	DLL *pDLL = new DLL;
	pDLL->pFileData = data;
	pDLL->dwFileSize = size;
	return pDLL;
}


bool CopySections(DLL* pDLL)
{
	int numberOfSections;
	PIMAGE_SECTION_HEADER pSections;
	DWORD sizeOfImage;
	if (pDLL->is64bit)
	{
		numberOfSections = pDLL->pNtHeaders64->FileHeader.NumberOfSections;
		pSections = IMAGE_FIRST_SECTION(pDLL->pNtHeaders64);
		sizeOfImage = pDLL->pNtHeaders64->OptionalHeader.SizeOfImage;
	}
	else
	{
		numberOfSections = pDLL->pNtHeaders32->FileHeader.NumberOfSections;
		pSections = IMAGE_FIRST_SECTION(pDLL->pNtHeaders32);
		sizeOfImage = pDLL->pNtHeaders32->OptionalHeader.SizeOfImage;
	}

	// calculate allocation requirements
	DWORD dwAlocationStart = 0xffffffff;
	DWORD dwAlocationEnd = 0;
	pDLL->Sections.reserve(numberOfSections);
	for (int sectionIdx = 0; sectionIdx < numberOfSections; sectionIdx++)
	{
		PIMAGE_SECTION_HEADER pSectionHeader = pSections + sectionIdx;
		if (pSectionHeader->Characteristics & IMAGE_SCN_MEM_DISCARDABLE)
			continue;
		DWORD dwSectionStart = pSectionHeader->VirtualAddress & ~(PAGE_SIZE-1);
		DWORD dwSectionEnd = (pSectionHeader->VirtualAddress + max(pSectionHeader->Misc.VirtualSize, pSectionHeader->SizeOfRawData) + PAGE_SIZE-1) & ~(PAGE_SIZE-1);
		if (dwAlocationEnd < dwSectionEnd)
			dwAlocationEnd = dwSectionEnd;
	}
	pDLL->dwAllocationSize = dwAlocationEnd;

	if (pDLL->dwAllocationSize / PAGE_SIZE > 65535)
	{
		fprintf(stderr, "error: DLL image is too big\n");
		return false;
	}

	// allocate image buffer
	pDLL->pImageBase = new BYTE[max(sizeOfImage, pDLL->dwAllocationSize)];
	if (pDLL->pImageBase == NULL)
	{
		fprintf(stderr, "error: could not allocate buffer for DLL image\n");
		return false;
	}
	memset(pDLL->pImageBase, 0, pDLL->dwAllocationSize);

	// copy sections
	for (int sectionIdx = 0; sectionIdx < numberOfSections; sectionIdx++)
	{
		PIMAGE_SECTION_HEADER pSectionHeader = pSections + sectionIdx;
		if (pSectionHeader->SizeOfRawData == 0)
			continue;
		if (pSectionHeader->VirtualAddress > sizeOfImage ||
			pSectionHeader->VirtualAddress + pSectionHeader->SizeOfRawData > sizeOfImage)
		{
			fprintf(stderr, "error: section data overflows image\n");
			return false;
		}
		memcpy(pDLL->pImageBase + pSectionHeader->VirtualAddress, pDLL->pFileData + pSectionHeader->PointerToRawData, pSectionHeader->SizeOfRawData);
	}

	return true;
}


bool BuildRelocations(DLL* pDLL)
{
	PIMAGE_DATA_DIRECTORY pRelocDir = 
		pDLL->is64bit ? &pDLL->pNtHeaders64->OptionalHeader.DataDirectory[IMAGE_DIRECTORY_ENTRY_BASERELOC]
					  :	&pDLL->pNtHeaders32->OptionalHeader.DataDirectory[IMAGE_DIRECTORY_ENTRY_BASERELOC];

	std::set<DWORD> relocations;

	// get relocations
	if (pRelocDir->Size > 0)
	{
		for (PIMAGE_BASE_RELOCATION pRelocation = (PIMAGE_BASE_RELOCATION)(pDLL->pImageBase + pRelocDir->VirtualAddress);
			 (DWORD)((PBYTE)pRelocation - (PBYTE)(pDLL->pImageBase + pRelocDir->VirtualAddress)) < pRelocDir->Size;
			 pRelocation = (PIMAGE_BASE_RELOCATION)(((PBYTE)pRelocation) + pRelocation->SizeOfBlock))
		{
			DWORD entriesCount = ((pRelocation->SizeOfBlock - sizeof(IMAGE_BASE_RELOCATION)) / 2);
			PWORD relInfo = (PWORD)((PBYTE)pRelocation + sizeof(IMAGE_BASE_RELOCATION));
			for (DWORD i = 0; i < entriesCount; i++, relInfo++)
			{
				int type, offset;

				// the upper 4 bits define the type of relocation
				type = *relInfo >> 12;

				// the lower 12 bits define the offset
				offset = *relInfo & 0xfff;
				
				switch (type)
				{
				case IMAGE_REL_BASED_ABSOLUTE:
					// skip relocation
					break;

				case IMAGE_REL_BASED_HIGHLOW:
					// change complete 32 bit address
					if (!pDLL->is64bit)
					{
						DWORD dwRelocationAddr = pRelocation->VirtualAddress + offset;
						std::set<DWORD>::iterator it = relocations.upper_bound(dwRelocationAddr);
						if (it != relocations.end())
						{
							if ((*it - dwRelocationAddr) < 4 || (it != relocations.begin() && (dwRelocationAddr - *(--it)) < 4))
							{
								fprintf(stderr, "error: found overlapping relocations at %08X\n", dwRelocationAddr);
								return false;
							}
						} 
						else if (!relocations.empty() && (dwRelocationAddr - *relocations.rbegin()) < 4)
						{
							fprintf(stderr, "error: found overlapping relocations at %08X\n", dwRelocationAddr);
							return false;
						}
						if (dwRelocationAddr > pDLL->dwAllocationSize - 4)
						{
							fprintf(stderr, "error: found out of range relocation\n");
							return false;
						}
						relocations.insert(dwRelocationAddr);
						*(DWORD*)(pDLL->pImageBase + dwRelocationAddr) -= 
							pDLL->is64bit ? (DWORD)pDLL->pNtHeaders64->OptionalHeader.ImageBase : pDLL->pNtHeaders32->OptionalHeader.ImageBase;
					}
					else
					{
						fprintf(stderr, "error: found unsupported IMAGE_REL_BASED_HIGHLOW relocation in 64-bit module\n");
						return false;
					}
					break;

				case IMAGE_REL_BASED_DIR64:
					// change complete 64 bit address
					if (pDLL->is64bit)
					{
						DWORD dwRelocationAddr = pRelocation->VirtualAddress + offset;
						std::set<DWORD>::iterator it = relocations.upper_bound(dwRelocationAddr);
						if (it != relocations.end())
						{
							if ((*it - dwRelocationAddr) < 8 || (it != relocations.begin() && (dwRelocationAddr - *(--it)) < 8))
							{
								fprintf(stderr, "error: found overlapping relocations at %08X\n", dwRelocationAddr);
								return false;
							}
						} 
						else if (!relocations.empty() && (dwRelocationAddr - *relocations.rbegin()) < 8)
						{
							fprintf(stderr, "error: found overlapping relocations at %08X\n", dwRelocationAddr);
							return false;
						}
						if (dwRelocationAddr > pDLL->dwAllocationSize - 8)
						{
							fprintf(stderr, "error: found out of range relocation\n");
							return false;
						}
						relocations.insert(dwRelocationAddr);
						*(DWORD64*)(pDLL->pImageBase + dwRelocationAddr) -= pDLL->pNtHeaders64->OptionalHeader.ImageBase;
					}
					else
					{
						fprintf(stderr, "error: found unsupported IMAGE_REL_BASED_DIR64 relocation in 32-bit module\n");
						return false;
					}
					break;

				
				default:
					fprintf(stderr, "error: found unsupported relocation type %d\n", type);
					return false;
				}
			}
		}
	}

	// split relocations to rellocation blocks
	std::vector<RELOCATION_BLOCK> relocBlocks;
	if (!relocations.empty())
	{
		relocBlocks.resize(1);
		RELOCATION_BLOCK *pRelocBlock = &relocBlocks[0];
		pRelocBlock->dwRelocationBase = (*relocations.begin() & ~0xff);
		pRelocBlock->dwRelocationEntries = 0;
		for (std::set<DWORD>::iterator it = relocations.begin(); it != relocations.end(); it++)
		{
			DWORD dwRelocation = *it;
			if ((dwRelocation & ~0xff) != pRelocBlock->dwRelocationBase)
			{
				relocBlocks.resize(relocBlocks.size()+1);
				pRelocBlock = &relocBlocks[relocBlocks.size()-1];
				pRelocBlock->dwRelocationBase = (dwRelocation & ~0xff);
				pRelocBlock->dwRelocationEntries = 0;
			}
			pRelocBlock->Relocations[pRelocBlock->dwRelocationEntries++] = dwRelocation & 0xff;
		}
	}

	// encode relocation blocks
	unsigned int lsfr1 = 0x6751DAFE, lsfr2 = GetTickCount();
	std::string relocStream;
	DWORD dwLastRelocationBase = 0;
	DWORD nBaseAndSize;
	for (std::vector<RELOCATION_BLOCK>::iterator it = relocBlocks.begin(); it != relocBlocks.end(); it++)
	{
		RELOCATION_BLOCK &relocBlock = *it;
		while (relocBlock.dwRelocationBase - dwLastRelocationBase >= 512 * 256) 
		{
			unsigned int nCnt = min(63, (unsigned int)(relocBlock.dwRelocationBase - dwLastRelocationBase)/(512*256));
			nBaseAndSize = (nCnt - 1) + 0xff80;
			relocStream.append((char*)&nBaseAndSize, 2);
			dwLastRelocationBase += nCnt * 512*256;
		}
		if (relocBlock.dwRelocationBase - dwLastRelocationBase >= 256)
		{
			nBaseAndSize = (relocBlock.dwRelocationEntries - 1) + ((relocBlock.dwRelocationBase - dwLastRelocationBase - 256) >> 1);
			relocStream.append((char*)&nBaseAndSize, 2);
		}
		else
		{
			nBaseAndSize = (BYTE)((relocBlock.dwRelocationEntries - 1) + ((relocBlock.dwRelocationBase - dwLastRelocationBase) >> 1)) + 0x40;
			relocStream.append((char*)&nBaseAndSize, 1);
		}
		for (unsigned int idx = 0; idx < relocBlock.dwRelocationEntries; idx++)
		{
			unsigned int swapIdx = (GetLSFRBits(6, lsfr1) ^ GetLSFRBits(6, lsfr2)) % relocBlock.dwRelocationEntries;
			std::swap(relocBlock.Relocations[idx], relocBlock.Relocations[swapIdx]);
		}
		relocStream.append((char*)relocBlock.Relocations, relocBlock.dwRelocationEntries);
		dwLastRelocationBase = relocBlock.dwRelocationBase + 256;
	}
	nBaseAndSize = 0xFFBF;
	relocStream.append((char*)&nBaseAndSize, 2);
	pDLL->Relocations = std::move(relocStream);

	return true;
}

bool CleanupRelocations(DLL* pDLL)
{
	PIMAGE_DATA_DIRECTORY pRelocDir = 
		pDLL->is64bit ? &pDLL->pNtHeaders64->OptionalHeader.DataDirectory[IMAGE_DIRECTORY_ENTRY_BASERELOC]
					  :	&pDLL->pNtHeaders32->OptionalHeader.DataDirectory[IMAGE_DIRECTORY_ENTRY_BASERELOC];

	// clean up relocations in image
	if (pRelocDir->Size > 0)
	{
		PIMAGE_BASE_RELOCATION pRelocation = (PIMAGE_BASE_RELOCATION)(pDLL->pImageBase + pRelocDir->VirtualAddress);
		while (pRelocation->VirtualAddress > 0)
		{
			PIMAGE_BASE_RELOCATION pNextRelocation = (PIMAGE_BASE_RELOCATION)(((PBYTE)pRelocation) + pRelocation->SizeOfBlock);
			memset(pRelocation, 0, pRelocation->SizeOfBlock);
			pRelocation = pNextRelocation;
		}
	}

	return true;
}


bool BuildImportsTable(DLL* pDLL, IMPORT_TRANSLATION_TABLE * pImportsXlat)
{
	bool ret = true;
	PIMAGE_DATA_DIRECTORY pImportsDir = 
		pDLL->is64bit ? &pDLL->pNtHeaders64->OptionalHeader.DataDirectory[IMAGE_DIRECTORY_ENTRY_IMPORT]
					  :	&pDLL->pNtHeaders32->OptionalHeader.DataDirectory[IMAGE_DIRECTORY_ENTRY_IMPORT];

	std::map<ULONGLONG, IMPORT_TABLE_ENTRY> imports;

	// build imports list
	if (pImportsDir->Size > 0)
	{
		if (pDLL->is64bit)
		{
			for (PIMAGE_IMPORT_DESCRIPTOR pImportDescriptor = (PIMAGE_IMPORT_DESCRIPTOR)(pDLL->pImageBase + pImportsDir->VirtualAddress);
				 !IsBadReadPtr(pImportDescriptor, sizeof(IMAGE_IMPORT_DESCRIPTOR)) && pImportDescriptor->Name;
				 pImportDescriptor++)
			{
				ULONGLONG nImportReference;
				ULONGLONG *pThunkReference;
				if (pImportDescriptor->OriginalFirstThunk)
					pThunkReference = (ULONGLONG*) (pDLL->pImageBase + pImportDescriptor->OriginalFirstThunk);
				else
					pThunkReference = (ULONGLONG*) (pDLL->pImageBase + pImportDescriptor->FirstThunk);
				for (nImportReference = pImportDescriptor->FirstThunk; 
					 *pThunkReference; 
					 pThunkReference++, nImportReference += sizeof(ULONGLONG))
				{
					IMPORT_TABLE_ENTRY import;
					if (IMAGE_SNAP_BY_ORDINAL64(*pThunkReference))
					{
						import.dwImportOrdinal = IMAGE_ORDINAL64(*pThunkReference);
						import.ImportName.clear();
					}
					else
					{
						PIMAGE_IMPORT_BY_NAME pImportByName = (PIMAGE_IMPORT_BY_NAME) (pDLL->pImageBase + (*pThunkReference));
						import.ImportName = (char*)&pImportByName->Name;
						import.dwImportOrdinal = (DWORD)-1;
					}
					import.ModuleName = (char*)(pDLL->pImageBase + pImportDescriptor->Name);
					import.SymbolName = import.ImportName;
					bool bDllImport = !import.ModuleName.empty() && import.ModuleName != ".exe.";
					if (pImportsXlat)
					{
						auto it = bDllImport ? 
							pImportsXlat->dllImports.find(import.ImportName) : pImportsXlat->imports.find(import.ImportName);
						if (it == (bDllImport ? pImportsXlat->dllImports.end() : pImportsXlat->imports.end()))
						{
							fprintf(stderr, "error: missing mangled name for%s %s\n", bDllImport ? " DLL import" : " import", import.ImportName.c_str());					
							ret = false;
						}
						else
						{
							import.SymbolName = it->second;
						}
					}
					imports[nImportReference] = import;
				}
			}
		}
		else
		{
			for (PIMAGE_IMPORT_DESCRIPTOR pImportDescriptor = (PIMAGE_IMPORT_DESCRIPTOR)(pDLL->pImageBase + pImportsDir->VirtualAddress);
				 !IsBadReadPtr(pImportDescriptor, sizeof(IMAGE_IMPORT_DESCRIPTOR)) && pImportDescriptor->Name;
				 pImportDescriptor++)
			{
				DWORD dwImportReference;
				DWORD *pThunkReference;
				if (pImportDescriptor->OriginalFirstThunk)
					pThunkReference = (DWORD*) (pDLL->pImageBase + pImportDescriptor->OriginalFirstThunk);
				else
					pThunkReference = (DWORD*) (pDLL->pImageBase + pImportDescriptor->FirstThunk);
				for (dwImportReference = pImportDescriptor->FirstThunk; 
					 *pThunkReference; 
					 pThunkReference++, dwImportReference += sizeof(DWORD))
				{
					IMPORT_TABLE_ENTRY import;
					if (IMAGE_SNAP_BY_ORDINAL32(*pThunkReference))
					{
						import.dwImportOrdinal = IMAGE_ORDINAL32(*pThunkReference);
						import.ImportName.clear();
					}
					else
					{
						PIMAGE_IMPORT_BY_NAME pImportByName = (PIMAGE_IMPORT_BY_NAME) (pDLL->pImageBase + (*pThunkReference));
						import.ImportName = (char*)&pImportByName->Name;
						import.dwImportOrdinal = (DWORD)-1;
					}
					import.ModuleName = (char*)(pDLL->pImageBase + pImportDescriptor->Name);
					import.SymbolName = import.ImportName;
					bool bDllImport = !import.ModuleName.empty() && import.ModuleName != ".exe.";
					if (pImportsXlat)
					{
						auto it = bDllImport ? 
							pImportsXlat->dllImports.find(import.ImportName) : pImportsXlat->imports.find(import.ImportName);
						if (it == (bDllImport ? pImportsXlat->dllImports.end() : pImportsXlat->imports.end()))
						{
							fprintf(stderr, "error: missing mangled name for%s %s\n", bDllImport ? " DLL import" : " import", import.ImportName.c_str());					
							ret = false;
						}
						else
						{
							import.SymbolName = it->second;
						}
					}
					imports[dwImportReference] = import;
				}
			}
		}
	}

	// clean up imports in image
	if (pImportsDir->Size > 0)
	{
		if (pDLL->is64bit)
		{
			for (PIMAGE_IMPORT_DESCRIPTOR pImportDescriptor = (PIMAGE_IMPORT_DESCRIPTOR)(pDLL->pImageBase + pImportsDir->VirtualAddress);
				 !IsBadReadPtr(pImportDescriptor, sizeof(IMAGE_IMPORT_DESCRIPTOR)) && pImportDescriptor->Name;
				 pImportDescriptor++)
			{
				char *pModuleName = (char*)(pDLL->pImageBase + pImportDescriptor->Name);
				ULONGLONG *pImportReference;
				ULONGLONG *pThunkReference;

				if (pImportDescriptor->OriginalFirstThunk)
					pThunkReference = (ULONGLONG*) (pDLL->pImageBase + pImportDescriptor->OriginalFirstThunk);
				else
					pThunkReference = (ULONGLONG*) (pDLL->pImageBase + pImportDescriptor->FirstThunk);
				for (pImportReference = (ULONGLONG*) (pDLL->pImageBase + pImportDescriptor->FirstThunk); *pThunkReference; pThunkReference++, pImportReference++)
				{
					if (!IMAGE_SNAP_BY_ORDINAL64(*pThunkReference))
					{
						PIMAGE_IMPORT_BY_NAME pImportByName = (PIMAGE_IMPORT_BY_NAME) (pDLL->pImageBase + (*pThunkReference));
						memset((char*)&pImportByName->Name, 0, strlen((char*)&pImportByName->Name));
						memset(pImportByName, 0, sizeof(IMAGE_IMPORT_BY_NAME));
					}
					*pThunkReference = 0;
					*pImportReference = 0;
				}

				memset(pModuleName, 0, strlen(pModuleName));
				memset(pImportDescriptor, 0, sizeof(IMAGE_IMPORT_DESCRIPTOR));
			}
		} 
		else
		{
			for (PIMAGE_IMPORT_DESCRIPTOR pImportDescriptor = (PIMAGE_IMPORT_DESCRIPTOR)(pDLL->pImageBase + pImportsDir->VirtualAddress);
				 !IsBadReadPtr(pImportDescriptor, sizeof(IMAGE_IMPORT_DESCRIPTOR)) && pImportDescriptor->Name;
				 pImportDescriptor++)
			{
				char *pModuleName = (char*)(pDLL->pImageBase + pImportDescriptor->Name);
				DWORD *pImportReference;
				DWORD *pThunkReference;

				if (pImportDescriptor->OriginalFirstThunk)
					pThunkReference = (DWORD*) (pDLL->pImageBase + pImportDescriptor->OriginalFirstThunk);
				else
					pThunkReference = (DWORD*) (pDLL->pImageBase + pImportDescriptor->FirstThunk);
				for (pImportReference = (DWORD*) (pDLL->pImageBase + pImportDescriptor->FirstThunk); *pThunkReference; pThunkReference++, pImportReference++)
				{
					if (!IMAGE_SNAP_BY_ORDINAL32(*pThunkReference))
					{
						PIMAGE_IMPORT_BY_NAME pImportByName = (PIMAGE_IMPORT_BY_NAME) (pDLL->pImageBase + (*pThunkReference));
						memset((char*)&pImportByName->Name, 0, strlen((char*)&pImportByName->Name));
						memset(pImportByName, 0, sizeof(IMAGE_IMPORT_BY_NAME));
					}
					*pThunkReference = 0;
					*pImportReference = 0;
				}

				memset(pModuleName, 0, strlen(pModuleName));
				memset(pImportDescriptor, 0, sizeof(IMAGE_IMPORT_DESCRIPTOR));
			}
		}
	}

	// split imports into blocks
	for (std::map<ULONGLONG, IMPORT_TABLE_ENTRY>::iterator it = imports.begin(); it != imports.end(); it++)
	{
		IMPORT_TABLE_ENTRY &import = it->second;
		std::list<IMPORT_BLOCK> &blocks = (import.ModuleName.empty() || import.ModuleName == ".exe.") ? pDLL->Imports : pDLL->DllImports;
		if (blocks.empty() || blocks.back().Imports.size() == 255 || it->first > blocks.back().dwImportReference + blocks.back().Imports.size() * sizeof(DWORD))
		{
			IMPORT_BLOCK block;
			block.dwImportReference = (DWORD)it->first;
			block.Imports.push_back(import);
			blocks.push_back(block);
		}
		else
		{
			blocks.back().Imports.push_back(import);
		}
	}

	return ret;
}


bool BuildExportsTable(DLL* pDLL)
{
	PIMAGE_DATA_DIRECTORY pExportsDir = 
		pDLL->is64bit ? &pDLL->pNtHeaders64->OptionalHeader.DataDirectory[IMAGE_DIRECTORY_ENTRY_EXPORT]
					  :	&pDLL->pNtHeaders32->OptionalHeader.DataDirectory[IMAGE_DIRECTORY_ENTRY_EXPORT];

	// build export table
	if (pExportsDir->Size > 0)
	{
		PIMAGE_EXPORT_DIRECTORY pExports = (PIMAGE_EXPORT_DIRECTORY) (pDLL->pImageBase + pExportsDir->VirtualAddress);

		// get adresses of functions
		pDLL->Exports.resize(pExports->NumberOfFunctions);
		for (DWORD idx = 0; idx < pExports->NumberOfFunctions; idx++)
		{
			pDLL->Exports[idx].dwAddressOfFunction = ((DWORD*)(pDLL->pImageBase + pExports->AddressOfFunctions))[idx];
			pDLL->Exports[idx].dwExportOrdinal = pExports->Base + idx;
		}

		// assign names to exported functions
		for (DWORD idx = 0; idx < pExports->NumberOfNames; idx++)
		{
			WORD wFunctionIdx = ((WORD*)(pDLL->pImageBase + pExports->AddressOfNameOrdinals))[idx];
			char* pszFunctionName = (char*)(pDLL->pImageBase + ((DWORD*)(pDLL->pImageBase + pExports->AddressOfNames))[idx]);
			pDLL->Exports[wFunctionIdx].FunctionNames.push_back(pszFunctionName);
		}
	}

	// cleanup exports in image
	if (pExportsDir->Size > 0)
	{
		PIMAGE_EXPORT_DIRECTORY pExports = (PIMAGE_EXPORT_DIRECTORY) (pDLL->pImageBase + pExportsDir->VirtualAddress);

		char* pszModuleName = (char*)(pDLL->pImageBase + pExports->Name);
		if (pExports->Name)
			memset(pszModuleName, 0, strlen(pszModuleName));

		memset(pDLL->pImageBase + pExports->AddressOfFunctions, 0, pExports->NumberOfFunctions * sizeof(DWORD));

		if (pExports->NumberOfNames)
		{
			for (DWORD idx = 0; idx < pExports->NumberOfNames; idx++)
			{
				char* pszFunctionName = (char*)(pDLL->pImageBase + ((DWORD*)(pDLL->pImageBase + pExports->AddressOfNames))[idx]);
				memset(pszFunctionName, 0, strlen(pszFunctionName));
			}
			memset(pDLL->pImageBase + pExports->AddressOfNames, 0, pExports->NumberOfNames * sizeof(DWORD));
			memset(pDLL->pImageBase + pExports->AddressOfNameOrdinals, 0, pExports->NumberOfNames * sizeof(WORD));
		}
	}

	return true;
}

DWORD GetSectionProtection(DWORD dwCharacteristics)
{
	if (dwCharacteristics & IMAGE_SCN_MEM_EXECUTE)
	{
		if (dwCharacteristics & IMAGE_SCN_MEM_WRITE)
		{
			if (dwCharacteristics & IMAGE_SCN_MEM_READ)
				return PAGE_EXECUTE_READWRITE;
			else
				return PAGE_EXECUTE_WRITECOPY;
		}
		else
		{
			if (dwCharacteristics & IMAGE_SCN_MEM_READ)
				return PAGE_EXECUTE_READ;
			else
				return PAGE_EXECUTE;
		}
	}
	else
	{
		if (dwCharacteristics & IMAGE_SCN_MEM_WRITE)
		{
			if (dwCharacteristics & IMAGE_SCN_MEM_READ)
				return PAGE_READWRITE;
			else
				return PAGE_WRITECOPY;
		}
		else
		{
			if (dwCharacteristics & IMAGE_SCN_MEM_READ)
				return PAGE_READONLY;
			else
				return PAGE_NOACCESS;
		}
	}
}

BYTE* DiBitsDecode(BYTE *pModuleCode, BYTE* pOutputBuffer)
{
    BYTE TempBuffer[256];
    int DiBits[16];
    BYTE bSuspendInfo;
	int nInputBits;
	int cInputBits;
	
    DiBits[0x1] = 0 | (1 << 8) | (2 << 16) | (3 << 24); // 1 - 00 01 10 11 - 0 1 2 3
    DiBits[0x2] = 0 | (2 << 8) | (1 << 16) | (3 << 24); // 2 - 00 10 01 11 - 0 2 1 3
    DiBits[0x3] = 0 | (3 << 8) | (1 << 16) | (2 << 24); // 3 - 00 11 01 10 - 0 3 1 2
    DiBits[0x4] = 1 | (0 << 8) | (2 << 16) | (3 << 24); // 4 - 01 00 10 11 - 1 0 2 3
    DiBits[0x6] = 1 | (2 << 8) | (0 << 16) | (3 << 24); // 6 - 01 10 00 11 - 1 2 0 3
    DiBits[0x7] = 1 | (3 << 8) | (0 << 16) | (2 << 24); // 7 - 01 11 00 10 - 1 3 0 2
    DiBits[0x8] = 2 | (0 << 8) | (1 << 16) | (3 << 24); // 8 - 10 00 01 11 - 2 0 1 3
    DiBits[0x9] = 2 | (1 << 8) | (0 << 16) | (3 << 24); // 9 - 10 01 00 11 - 2 1 0 3
    DiBits[0xb] = 2 | (3 << 8) | (0 << 16) | (1 << 24); // B - 10 11 00 01 - 2 3 0 1
    DiBits[0xc] = 3 | (0 << 8) | (1 << 16) | (2 << 24); // C - 11 00 01 10 - 3 0 1 2
    DiBits[0xd] = 3 | (1 << 8) | (0 << 16) | (2 << 24); // D - 11 01 00 10 - 3 1 0 2
    DiBits[0xe] = 3 | (2 << 8) | (0 << 16) | (1 << 24); // E - 11 10 00 01 - 3 2 0 1

	// decode section data
	bSuspendInfo = 128;
	nInputBits = 0;
	cInputBits = 0;
	for (;;)
	{
		int cDiBitsBuffer;
		int cDiBits0, cDiBits1, cDiBits2, cDiBits3;
		BYTE *pDiBitsBuffer;
		DWORD nBytesCount;
    	BYTE *pRLEBuffer;

		// get bi-bits encoding and output bytes count
		if (nInputBits < 8)
		{
			cInputBits |= ((*pModuleCode++) << nInputBits);
			nInputBits += 8;
		}

		// decode output bytes count:
		// 0x0*-0x9* = 1-10 bytes
		// 0xA* = 16 bytes
		// 0xB* = 32 bytes
		// 0xC* = 64 bytes
		// 0xD* = 128 bytes
		// 0xE* = 192 bytes
		// 0xF* = 256 bytes
		// 0xFF = end-of-stream
		if ((cInputBits & 0xff) >= 0xf0)
		{
        	if (cInputBits == 0xff) 
        		break; // end of stream
        	nBytesCount = 0x100;
		}
		else if ((cInputBits & 0xff) >= 0xe0)
		{
        	nBytesCount = 0xc0;
		}
		else if ((cInputBits & 0xff) >= 0xa0)
		{
        	nBytesCount = 0x10 << (((cInputBits & 0xff) - 0xa0) >> 4);
		}
		else
		{
        	nBytesCount = ((cInputBits & 0xff) >> 4) + 1;
		}

		// setup di-bits
		cDiBits0 = DiBits[cInputBits & 0xf];
		cDiBits1 = (cDiBits0 >> 8) & 0xff;
		cDiBits2 = (cDiBits0 >> 16) & 0xff;
		cDiBits3 = (cDiBits0 >> 24) & 0xff;
		cDiBits0 &= 0xff;

		// consume 8 bits
		cInputBits >>= 8;
		nInputBits -= 8;

		// decode one di-bits block
		cDiBitsBuffer = 4;
		pDiBitsBuffer = TempBuffer + sizeof(TempBuffer) - nBytesCount;
		for (;;)
		{
			if (nInputBits < 3)
			{
				cInputBits |= ((*pModuleCode++) << nInputBits);
				nInputBits += 8;
			}

			if (cInputBits & 1)
			{
            	// 1 - cDiBits0
            	cDiBitsBuffer += cDiBits0;
            	cInputBits >>= 1;
				nInputBits--;
			}
			else if (!(cInputBits & 2))
			{
            	// 00 - cDiBits1
            	cDiBitsBuffer += cDiBits1;
            	cInputBits >>= 2;
				nInputBits -= 2;
			}
			else if (!(cInputBits & 4))
			{        
            	// 010 - cDiBits2
            	cDiBitsBuffer += cDiBits2;
            	cInputBits >>= 3;
				nInputBits -= 3;
			}
			else
			{
            	// 011 - cDiBits3
            	cDiBitsBuffer += cDiBits3;
            	cInputBits >>= 3;
				nInputBits -= 3;
			}

			if (cDiBitsBuffer > 255)
			{
				*pDiBitsBuffer++ = (BYTE)cDiBitsBuffer;
            	if (pDiBitsBuffer == TempBuffer + sizeof(TempBuffer))
            		break;
            	cDiBitsBuffer = 4;
			}
			else
			{
            	cDiBitsBuffer <<= 2;
			}
		}

		// setup RLE buffer
    	pRLEBuffer = TempBuffer + sizeof(TempBuffer) - nBytesCount;

       	// check if we suspended last time, clean up
       	if (bSuspendInfo > 128)
       	{
       		// suspended replicate run
       		DWORD nCount = (BYTE)(257 - bSuspendInfo);

       		// expand replicate run
       		DWORD dwReplicate = *pRLEBuffer++;
       		--nBytesCount;
			if (nCount & 1)
				*pOutputBuffer++ = (BYTE)dwReplicate;
			dwReplicate |= dwReplicate << 8;
			if (nCount & 2)
				*((WORD*&)pOutputBuffer)++ = (WORD)dwReplicate;
			dwReplicate |= dwReplicate << 16;
			nCount >>= 2;
			while (nCount--)
				*((DWORD*&)pOutputBuffer)++ = dwReplicate;

       		// not suspended now
       		bSuspendInfo = 128;
       	}
       	else if (bSuspendInfo < 128)
       	{
       		// suspended literal run
       		DWORD nCount = (BYTE)((bSuspendInfo & 0x7f) + 1);

       		if (nBytesCount < nCount)
       		{
       			// suspending again
       			bSuspendInfo = (BYTE)(nCount - nBytesCount - 1);
       			nCount = (BYTE)(nBytesCount);
       		}
       		else
       		{
       			// not suspended now
       			bSuspendInfo = 128;
       		}

       		// copy literals
       		nBytesCount -= nCount;
			if (nCount & 1)
				*pOutputBuffer++ = *pRLEBuffer++;
			if (nCount & 2)
				*((WORD*&)pOutputBuffer)++ = *((WORD*&)pRLEBuffer)++;
			nCount >>= 2;
			while (nCount--)
				*((DWORD*&)pOutputBuffer)++ = *((DWORD*&)pRLEBuffer)++;
       	}

       	// continue with ordinary RLE decompression
       	while (nBytesCount)
       	{
       		BYTE b = *pRLEBuffer++;
       		--nBytesCount;

       		if (b & 0x80)
       		{
       			// replicate run
       			if (nBytesCount)
       			{
       				DWORD nCount = (BYTE)(257 - b);

       				// expand replicate run
       				DWORD dwReplicate = *pRLEBuffer++;
       				--nBytesCount;
					if (nCount & 1)
						*pOutputBuffer++ = (BYTE)dwReplicate;
					dwReplicate |= dwReplicate << 8;
					if (nCount & 2)
						*((WORD*&)pOutputBuffer)++ = (WORD)dwReplicate;
					dwReplicate |= dwReplicate << 16;
					nCount >>= 2;
					while (nCount--)
						*((DWORD*&)pOutputBuffer)++ = dwReplicate;
       			}
       			else
       			{
       				// missing second byte of replicate run, prepare suspend
       				bSuspendInfo = b;
       			}
       		}
       		else
       		{
       			// literal run
       			DWORD nCount = (BYTE)((b & 0x7f) + 1);
       			if (nBytesCount < nCount)
       			{
       				// not enough literals, prepare suspend
       				bSuspendInfo = (BYTE)(nCount - nBytesCount - 1);
       				nCount = (BYTE)(nBytesCount);
       			}

       			// copy literals
       			nBytesCount -= nCount;
				if (nCount & 1)
					*pOutputBuffer++ = *pRLEBuffer++;
				if (nCount & 2)
					*((WORD*&)pOutputBuffer)++ = *((WORD*&)pRLEBuffer)++;
				nCount >>= 2;
				while (nCount--)
					*((DWORD*&)pOutputBuffer)++ = *((DWORD*&)pRLEBuffer)++;
       		}
       	}
	}

	return pOutputBuffer;
}

bool FinalizeSections(DLL* pDLL)
{
	DWORD nSectionsCount = 0;
	int numberOfSections;
	PIMAGE_SECTION_HEADER pSections;
	DWORD sizeOfImage;
	if (pDLL->is64bit)
	{
		numberOfSections = pDLL->pNtHeaders64->FileHeader.NumberOfSections;
		pSections = IMAGE_FIRST_SECTION(pDLL->pNtHeaders64);
		sizeOfImage = pDLL->pNtHeaders64->OptionalHeader.SizeOfImage;
	}
	else
	{
		numberOfSections = pDLL->pNtHeaders32->FileHeader.NumberOfSections;
		pSections = IMAGE_FIRST_SECTION(pDLL->pNtHeaders32);
		sizeOfImage = pDLL->pNtHeaders32->OptionalHeader.SizeOfImage;
	}

	// build sections table
	pDLL->Sections.reserve(numberOfSections);
	for (int sectionIdx = 0; sectionIdx < numberOfSections; sectionIdx++)
	{
		PIMAGE_SECTION_HEADER pSectionHeader = pSections + sectionIdx;
		if (pSectionHeader->Characteristics & IMAGE_SCN_MEM_DISCARDABLE)
			continue;

		pDLL->Sections.resize(++nSectionsCount);
		SECTION_TABLE_ENTRY &section = pDLL->Sections.back();

		if (opt_SingleSection)
		{
			section.dwSectionBase = pSectionHeader->VirtualAddress;
			section.dwSectionSize = pSectionHeader->Misc.VirtualSize;
			section.dwDataSize = pSectionHeader->Misc.VirtualSize;
		}
		else
		{
			section.dwSectionBase = pSectionHeader->VirtualAddress & ~(PAGE_SIZE-1);
			section.dwSectionSize = (max(pSectionHeader->Misc.VirtualSize, pSectionHeader->SizeOfRawData) + PAGE_SIZE-1) & ~(PAGE_SIZE-1);
			section.dwDataSize = (pSectionHeader->VirtualAddress - section.dwSectionBase) + pSectionHeader->SizeOfRawData;
		}
		section.dwProtection = GetSectionProtection(pSectionHeader->Characteristics);

		if (opt_LZMA)
		{
			CLzmaEncProps props;
			SRes res;
			Byte *dest;
			SizeT destLen;
			const Byte *src;
			SizeT srcLen;
			Byte propsEncoded[5];
			SizeT propsSize;
			int writeEndMark;
			ELzmaStatus status;

			// initialize LZMA properties
			LzmaEncProps_Init(&props);
			props.level = 9;
			LzmaEncProps_Normalize(&props);
			props.lc = 1;
			props.lp = 0;
			props.pb = 1;
			props.writeEndMark = 1;

			// compress data
			section.Data.resize(2 * section.dwDataSize + 2);
			dest = (Byte *)section.Data.data();
			destLen = section.Data.size();
			src = (const Byte*)(pDLL->pImageBase + section.dwSectionBase);
			srcLen = section.dwDataSize;
			propsSize = sizeof(propsEncoded);
			writeEndMark = 1;
			res = LzmaEncode(dest, &destLen, src, srcLen, &props, propsEncoded, &propsSize, writeEndMark, NULL, &g_Alloc, &g_Alloc);
			if (res != SZ_OK)
			{
				fprintf(stderr, "error: LZMA compression failed\n");
				return false;
			}
			section.Data[destLen] = propsEncoded[0];
			DWORD dicSize = propsEncoded[1] | ((UInt32)propsEncoded[2] << 8) | ((UInt32)propsEncoded[3] << 16) | ((UInt32)propsEncoded[4] << 24);
			section.Data[destLen+1] = 0;
			while ((DWORD)(1 << (section.Data[destLen+1])) < dicSize)
				section.Data[destLen+1]++;
			section.Data.resize(destLen + 2);

			// verify compressed data
			BYTE *pDecodedData = new BYTE[section.dwDataSize];
			memset(pDecodedData, 0, section.dwDataSize);
			dest = pDecodedData;
			destLen = section.dwDataSize;
			src = (const Byte*)section.Data.data();
			srcLen = section.Data.size();
			if (SZ_OK != LzmaDecode(pDecodedData, &destLen, src, &srcLen, propsEncoded, propsSize, LZMA_FINISH_END, &status, &g_Alloc) ||
				status != LZMA_STATUS_FINISHED_WITH_MARK ||
				memcmp(pDLL->pImageBase + section.dwSectionBase, pDecodedData, section.dwDataSize) != 0)
			{
				delete[] pDecodedData;
				fprintf(stderr, "error: LZMA compression failed\n");
				return false;
			}
			delete[] pDecodedData;
		}
		else if (opt_DiBits)
		{
			// RLE encode section data
			std::string rleData;
			RLEEncoder rle;
			rle.Encode(pDLL->pImageBase + section.dwSectionBase, section.dwDataSize);
			rleData.resize(rle.GetOutputSize());
			rle.WriteOutput((BYTE*)rleData.c_str());

			// di-bits encode section data
			DiBitsEncoder dibits;
			dibits.Encode((BYTE*)rleData.c_str(), rleData.size());
			section.Data.resize(dibits.GetOutputSize());
			dibits.WriteOutput((BYTE*)section.Data.c_str());

			// verify encoded data
			BYTE *pDecodedData = new BYTE[section.dwDataSize];
			memset(pDecodedData, 0, section.dwDataSize);
			if (DiBitsDecode((BYTE*)section.Data.c_str(), pDecodedData) != pDecodedData + section.dwDataSize ||
				memcmp(pDLL->pImageBase + section.dwSectionBase, pDecodedData, section.dwDataSize) != 0)
			{
				delete[] pDecodedData;
				fprintf(stderr, "error: Di-Bits compression failed\n");
				return false;
			}
			delete[] pDecodedData;
		}
		else
		{
			section.Data.assign((const char*)(pDLL->pImageBase + section.dwSectionBase), section.dwDataSize);
		}
	}

	return true;
}


bool ParseDll(DLL* pDLL, IMPORT_TRANSLATION_TABLE * pImportsXlat)
{
	PIMAGE_DOS_HEADER pDosHeader = NULL;
	PIMAGE_NT_HEADERS32 pNtHeaders32 = NULL;
	PIMAGE_NT_HEADERS64 pNtHeaders64 = NULL;

	// check file headers
	pDosHeader = (PIMAGE_DOS_HEADER)pDLL->pFileData;
	if (pDosHeader->e_magic != IMAGE_DOS_SIGNATURE)
	{
		fprintf(stderr, "error: invalid DLL file format\n");
		return false;
	}
	pDLL->pDosHeader = pDosHeader;

	PIMAGE_NT_HEADERS pHdr = ((DWORD)pDosHeader->e_lfanew <= pDLL->dwFileSize - sizeof(IMAGE_NT_HEADERS)) ? (PIMAGE_NT_HEADERS)((PBYTE)pDLL->pFileData + pDosHeader->e_lfanew) : NULL;
	if (!pHdr || pHdr->Signature != IMAGE_NT_SIGNATURE)
	{
		fprintf(stderr, "error: no PE header found in DLL\n");
		return false;
	}
	pDLL->is64bit = (pHdr->FileHeader.Machine == IMAGE_FILE_MACHINE_AMD64);
	if (pDLL->is64bit)
	{
		opt_64bit = true;
		pDLL->pNtHeaders32 = NULL;
		pDLL->pNtHeaders64 = (PIMAGE_NT_HEADERS64)pHdr;
	}
	else
	{
		pDLL->pNtHeaders32 = (PIMAGE_NT_HEADERS32)pHdr;
		pDLL->pNtHeaders64 = NULL;
	}

	// check section alignment
	if(pDLL->is64bit && pDLL->pNtHeaders64->OptionalHeader.SectionAlignment != 4096)
	{
		fprintf(stderr, "error: unsupported section alignment (%d bytes), sections must be aligned to 4096 bytes\n", pDLL->pNtHeaders64->OptionalHeader.SectionAlignment);
		return false;
	}
	else if(!pDLL->is64bit && pDLL->pNtHeaders32->OptionalHeader.SectionAlignment != 4096)
	{
		fprintf(stderr, "error: unsupported section alignment (%d bytes), sections must be aligned to 4096 bytes\n", pDLL->pNtHeaders32->OptionalHeader.SectionAlignment);
		return false;
	}

	// implicit TLS are not supported
	if (pDLL->is64bit && pDLL->pNtHeaders64->OptionalHeader.DataDirectory[IMAGE_DIRECTORY_ENTRY_TLS].Size != 0)
	{
		fprintf(stderr, "error: DLL contains implicit TLS\n");
		return false;
	}
	else if (!pDLL->is64bit && pDLL->pNtHeaders32->OptionalHeader.DataDirectory[IMAGE_DIRECTORY_ENTRY_TLS].Size != 0)
	{
		fprintf(stderr, "error: DLL contains implicit TLS\n");
		return false;
	}

	if (!CopySections(pDLL))
		return false;

	if (opt_SingleSection)
	{
		if (pDLL->is64bit)
		{
			if (pDLL->pNtHeaders64->FileHeader.NumberOfSections == 2)
			{
				for (int sectionIdx = 0; sectionIdx < pDLL->pNtHeaders64->FileHeader.NumberOfSections; sectionIdx++)
				{
					PIMAGE_SECTION_HEADER pSectionHeader = IMAGE_FIRST_SECTION(pDLL->pNtHeaders64) + sectionIdx;
					std::string sectionName;
					sectionName = std::string((const char*)pSectionHeader->Name, strnlen((const char*)pSectionHeader->Name, 8));
					if (sectionName != ".text" && sectionName != ".pdata")
					{
						fprintf(stderr, "error: DLL must contain only .text and .pdata sections\n");
						return false;
					}
				}
			}
			else if (pDLL->pNtHeaders64->FileHeader.NumberOfSections > 1)
			{
				fprintf(stderr, "error: DLL must contain only one section\n");
				return false;
			}
			if (pDLL->pNtHeaders64->OptionalHeader.DataDirectory[IMAGE_DIRECTORY_ENTRY_BASERELOC].Size > 0)
			{
				fprintf(stderr, "error: DLL must not contain relocations\n");
				return false;
			}
			if (pDLL->pNtHeaders64->OptionalHeader.DataDirectory[IMAGE_DIRECTORY_ENTRY_IMPORT].Size > 0)
			{
				fprintf(stderr, "error: DLL must not contain imports\n");
				return false;
			}
			if (pDLL->pNtHeaders64->OptionalHeader.DataDirectory[IMAGE_DIRECTORY_ENTRY_EXPORT].Size > 0)
			{
				fprintf(stderr, "error: DLL must not contain exports\n");
				return false;
			}
			if (pDLL->pNtHeaders64->OptionalHeader.AddressOfEntryPoint != 0 &&
				pDLL->pNtHeaders64->OptionalHeader.AddressOfEntryPoint != ((PIMAGE_SECTION_HEADER)IMAGE_FIRST_SECTION(pDLL->pNtHeaders64))->VirtualAddress)
			{
				fprintf(stderr, "error: DLL entry point must point to start of section\n");
				return false;
			}
		}
		else
		{
			if (pDLL->pNtHeaders32->FileHeader.NumberOfSections == 2)
			{
				for (int sectionIdx = 0; sectionIdx < pDLL->pNtHeaders32->FileHeader.NumberOfSections; sectionIdx++)
				{
					PIMAGE_SECTION_HEADER pSectionHeader = IMAGE_FIRST_SECTION(pDLL->pNtHeaders32) + sectionIdx;
					std::string sectionName;
					sectionName = std::string((const char*)pSectionHeader->Name, strnlen((const char*)pSectionHeader->Name, 8));
					if (sectionName != ".text" && sectionName != ".rdata")
					{
						fprintf(stderr, "error: DLL must contain only .text and .rdata sections\n");
						return false;
					}
				}
			}
			else if (pDLL->pNtHeaders32->FileHeader.NumberOfSections > 1)
			{
				fprintf(stderr, "error: DLL must contain only one section\n");
				return false;
			}
			if (pDLL->pNtHeaders32->OptionalHeader.DataDirectory[IMAGE_DIRECTORY_ENTRY_BASERELOC].Size > 0)
			{
				fprintf(stderr, "error: DLL must not contain relocations\n");
				return false;
			}
			if (pDLL->pNtHeaders32->OptionalHeader.DataDirectory[IMAGE_DIRECTORY_ENTRY_IMPORT].Size > 0)
			{
				fprintf(stderr, "error: DLL must not contain imports\n");
				return false;
			}
			if (pDLL->pNtHeaders32->OptionalHeader.DataDirectory[IMAGE_DIRECTORY_ENTRY_EXPORT].Size > 0)
			{
				fprintf(stderr, "error: DLL must not contain exports\n");
				return false;
			}
			if (pDLL->pNtHeaders32->OptionalHeader.AddressOfEntryPoint != 0 &&
				pDLL->pNtHeaders32->OptionalHeader.AddressOfEntryPoint != ((PIMAGE_SECTION_HEADER)IMAGE_FIRST_SECTION(pDLL->pNtHeaders32))->VirtualAddress)
			{
				fprintf(stderr, "error: DLL entry point must point to start of section\n");
				return false;
			}
		}
	}
	else
	{
		if (!BuildRelocations(pDLL))
			return false;
		if (!CleanupRelocations(pDLL))
			return false;

		if (!BuildImportsTable(pDLL, pImportsXlat))
			return false;

		if (!BuildExportsTable(pDLL))
			return false;
	}

	if (!FinalizeSections(pDLL))
		return false;

	return true;
}


COFF* LoadCOFFToMemory(const char* pszObjFilePath)
{
	int fd;
	size_t size;
	PBYTE data = NULL;
	
	if (_sopen_s(&fd, pszObjFilePath, _O_BINARY | _O_RDONLY, _SH_DENYNO, 0))
	{
		fprintf(stderr, "error: could not open COFF file \"%s\"\n", pszObjFilePath);
		return NULL;
	}

	size = _filelength(fd);
	data = new BYTE[size];
	if (data == NULL)
	{
		_close(fd);
		fprintf(stderr, "error: could not allocate buffer for file \"%s\"\n", pszObjFilePath);
		return NULL;
	}
	if (_read(fd, data, size) != size)
	{
		delete data;
		_close(fd);
		fprintf(stderr, "error: could not read file \"%s\"\n", pszObjFilePath);
		return NULL;
	}
	_close(fd);

	COFF *pCOFF = new COFF;
	pCOFF->pFileData = data;
	pCOFF->dwFileSize = size;
	return pCOFF;
}


bool ParseCOFF(COFF* pCOFF)
{
	PIMAGE_FILE_HEADER pCoffHeader;

	// check file header
	if (pCOFF->dwFileSize < sizeof(IMAGE_FILE_HEADER))
	{
		fprintf(stderr, "error: invalid format of COFF file\n");
		return false;
	}
	pCoffHeader = (PIMAGE_FILE_HEADER)pCOFF->pFileData;
	if (pCoffHeader->Machine != IMAGE_FILE_MACHINE_UNKNOWN &&
		pCoffHeader->Machine != IMAGE_FILE_MACHINE_I386 &&
		pCoffHeader->Machine != IMAGE_FILE_MACHINE_AMD64)
	{
		fprintf(stderr, "error: unupported machine type (%04X) in COFF file\n", (int)pCoffHeader->Machine);
		return false;
	}

	pCOFF->is64bit = (pCoffHeader->Machine == IMAGE_FILE_MACHINE_AMD64);

	// set string table pointer
	const char *pStringTable = (const char*)((PIMAGE_SYMBOL)((ULONG_PTR)pCoffHeader + pCoffHeader->PointerToSymbolTable) + pCoffHeader->NumberOfSymbols);

	// load symbols
	pCOFF->Symbols.resize(pCoffHeader->NumberOfSymbols);
	for (DWORD symbolIdx = 0; symbolIdx < pCoffHeader->NumberOfSymbols; symbolIdx++)
	{
		PIMAGE_SYMBOL pSymbolHeader = (PIMAGE_SYMBOL)((ULONG_PTR)pCoffHeader + pCoffHeader->PointerToSymbolTable) + symbolIdx;
		symbolIdx += pSymbolHeader->NumberOfAuxSymbols;

		COFF_SYMBOL &symbol = pCOFF->Symbols[symbolIdx];

		if (pSymbolHeader->N.Name.Short == 0)
			symbol.name = pStringTable + pSymbolHeader->N.Name.Long;
		else
			symbol.name = std::string((const char*)pSymbolHeader->N.ShortName, strnlen((const char*)pSymbolHeader->N.ShortName, 8));
		symbol.value = pSymbolHeader->Value;
		symbol.sectionNumber = pSymbolHeader->SectionNumber;
		symbol.type = pSymbolHeader->Type;
		symbol.storageClass = pSymbolHeader->StorageClass;
		symbol.auxDataSize = pSymbolHeader->NumberOfAuxSymbols * sizeof(IMAGE_SYMBOL);
		symbol.auxData = pSymbolHeader->NumberOfAuxSymbols ? (PIMAGE_AUX_SYMBOL)(pSymbolHeader + 1) : nullptr;
	}

	// load sections
	pCOFF->Sections.resize(pCoffHeader->NumberOfSections);
	for (WORD sectionIdx = 0; sectionIdx < pCoffHeader->NumberOfSections; sectionIdx++)
	{
		COFF_SECTION &section = pCOFF->Sections[sectionIdx];
		PIMAGE_SECTION_HEADER pSectionHeader = 	(PIMAGE_SECTION_HEADER)((BYTE*)pCoffHeader +
			sizeof(IMAGE_FILE_HEADER) + pCoffHeader->SizeOfOptionalHeader) + sectionIdx;

		if (pSectionHeader->Name[0] == '/')
			section.name = pStringTable + atoi((const char*)pSectionHeader->Name + 1);
		else
			section.name = std::string((const char*)pSectionHeader->Name, strnlen((const char*)pSectionHeader->Name, 8));
		section.virtualSize = pSectionHeader->Misc.VirtualSize;
		section.virtualAddress = pSectionHeader->VirtualAddress;
		section.rawDataSize = pSectionHeader->SizeOfRawData;
		section.rawData = (BYTE*)pCoffHeader + pSectionHeader->PointerToRawData;
		section.relocationsCount = pSectionHeader->NumberOfRelocations;
		section.relocations = pSectionHeader->NumberOfRelocations ? 
			(PIMAGE_RELOCATION)((BYTE*)pCoffHeader + pSectionHeader->PointerToRelocations) : nullptr;
		section.characteristics = pSectionHeader->Characteristics;
	}

	return true;
}


bool BuildImportsTranslationTable(COFF* pCOFF, IMPORT_TRANSLATION_TABLE * pXlat)
{
	for (COFF_SECTION & section : pCOFF->Sections)
	{
		bool bDllImport = false;
		const char *svmSymbol;
		const char *mangledName;

		if (strncmp(section.name.c_str(), ".isxt$", 6))
			continue;

		svmSymbol = section.name.c_str() + 6;
		if (!strncmp(svmSymbol, "_dll_", 5))
		{
			svmSymbol += 5;
			bDllImport = true;
		}

		if (section.relocationsCount == 0)
		{
			fprintf(stderr, "error: missing mangled name for symbol %s\n", svmSymbol);
			return false;
		}
		if (section.relocationsCount != 1)
		{
			fprintf(stderr, "error: multiple mangled names defined for symbol %s\n", svmSymbol);
			return false;
		}

		if (section.relocations[0].SymbolTableIndex >= pCOFF->Symbols.size())
		{
			fprintf(stderr, "error: missing mangled name for symbol %s\n", svmSymbol);
			return false;
		}
		mangledName = pCOFF->Symbols[section.relocations[0].SymbolTableIndex].name.c_str();

		auto & table = bDllImport ? pXlat->dllImports : pXlat->imports;
		if (table.find(svmSymbol) != table.end())
		{
			fprintf(stderr, "error: multiple mangled names defined for symbol %s\n", svmSymbol);
			return false;
		}
		table[svmSymbol] = mangledName;
	}

	return true;
}


std::string OutputMemberName()
{
	static unsigned int memberIdx = 0;
	char buf[32];
	sprintf_s<_countof(buf)>(buf, "_%u", ++memberIdx);
	return buf;
}


std::string OutputBYTE(BYTE num)
{
	char buf[32];
	if (opt_Asm)
		sprintf_s<_countof(buf)>(buf, "0%02XH", (int)num);
	else
		sprintf_s<_countof(buf)>(buf, "0x%02x", (int)num);
	return buf;
}

std::string OutputWORD(WORD num)
{
	char buf[32];
	if (opt_Asm)
		sprintf_s<_countof(buf)>(buf, "0%04XH", (int)num);
	else
		sprintf_s<_countof(buf)>(buf, "0x%04x", (int)num);
	return buf;
}

std::string OutputDWORD(DWORD num)
{
	char buf[32];
	if (opt_Asm)
		sprintf_s<_countof(buf)>(buf, "0%08XH", (int)num);
	else
		sprintf_s<_countof(buf)>(buf, "0x%08x", (int)num);
	return buf;
}

std::string OutputQWORD(DWORD64 num)
{
	char buf[32];
	if (opt_Asm)
		sprintf_s<_countof(buf)>(buf, "0%016llXH", (int64_t)num);
	else
		sprintf_s<_countof(buf)>(buf, "0x%016llx", (int64_t)num);
	return buf;
}

std::string OutputOFFSET(uint64_t num)
{
	char buf[32];
	if (opt_64bit)
	{
		if (opt_Asm)
			sprintf_s<_countof(buf)>(buf, "0%016llXH", (long long)num);
		else
			sprintf_s<_countof(buf)>(buf, "0x%016llx", (long long)num);
	}
	else
	{
		if (opt_Asm)
			sprintf_s<_countof(buf)>(buf, "0%08XH", (int)num);
		else
			sprintf_s<_countof(buf)>(buf, "0x%08x", (int)num);
	}
	return buf;
}

std::string OutputNumber(size_t num)
{
	char buf[32];
	sprintf_s<_countof(buf)>(buf, "%u", (int)num);
	return buf;
}

std::string OutputData(int tabs, const BYTE* data, size_t size)
{
	char buf[256];
	std::string str;
	while (size)
	{
		str.append(tabs, '\t');
		if (size >= 16)
		{
			if (opt_Asm)
			{
				const DWORD *d = (const DWORD *)data;
				sprintf_s<_countof(buf)>(buf, "DD 0%08XH,0%08XH,0%08XH,0%08XH", (int)d[0], (int)d[1], (int)d[2], (int)d[3]);
			}
			else
			{
				sprintf_s<_countof(buf)>(buf, "0x%02x,0x%02x,0x%02x,0x%02x,0x%02x,0x%02x,0x%02x,0x%02x,0x%02x,0x%02x,0x%02x,0x%02x,0x%02x,0x%02x,0x%02x,0x%02x", 
					(int)data[0], (int)data[1], (int)data[2], (int)data[3], (int)data[4], (int)data[5], (int)data[6], (int)data[7],
					(int)data[8], (int)data[9], (int)data[10], (int)data[11], (int)data[12], (int)data[13], (int)data[14], (int)data[15]);
			}
			str += buf;
			data += 16;
			size -= 16;
		}
		else
		{
			if (opt_Asm)
				str += "DB ";
			while (size)
			{
				if (opt_Asm)
					sprintf_s<_countof(buf)>(buf, "0%02XH", (int)*data);
				else
					sprintf_s<_countof(buf)>(buf, "0x%02x", (int)*data);
				str += buf;
				++data;
				if (--size)
					str += ",";
			}
		}
		if (size && !opt_Asm)
			str += ",\n";
		else
			str += "\n";
	}
	return str;
}

uint64_t GenOFFSET(unsigned int & lcg)
{
	if (opt_64bit)
	{
		uint64_t offset;
		offset = ((uint64_t)(LCG_69069_1(lcg) >> 2) << 32ULL);
		offset |= LCG_69069_1(lcg);
		offset |= 0x8000000000000000ULL;
		return offset;
	}
	else
	{
		return ((LCG_69069_1(lcg) >> 2) | 0x80000000);
	}
}
#define SYMBOL_OFFSET OutputOFFSET(GenOFFSET(lcg))


bool ExportDLL(
	DLL *pDLL,
	const char* pszModuleHeaderPath,
	const char* pszModuleDataPath,
	const char* pszSymbolName)
{
	static const struct {
		int			v;
		const char* i;
	} opcodes[15] = { 
		{ 0xb8, "mov eax" },
		{ 0xb9, "mov ecx" },
		{ 0xba, "mov edx" },
		{ 0xbb, "mov ebx" },
		{ 0xbc, "mov esp" },
		{ 0xbd, "mov ebp" },
		{ 0xbe, "mov esi" },
		{ 0xbf, "mov edi" },
		{ 0x05, "add eax" },
		{ 0x0d, "or eax"  },
		{ 0x15, "adc eax" },
		{ 0x1d, "sbb eax" },
		{ 0x25, "and eax" },
		{ 0x2d, "sub eax" },
		{ 0x35, "xor eax" } 
		// 0x3d - "cmp eax" - used as delimiter
	};
	int fd;
	std::string header;
	std::string decl;
	std::string data;
	bool ret = true;

	header += "#pragma once\n\n";
	header += "extern \"C\" void " + std::string(pszSymbolName) + "_DATA();\n\n";
	header += "#pragma pack(push,1)\nstruct " + std::string(pszSymbolName) + " {\n";
	header += "\tstatic const void*" + (opt_64bit ? std::string() : std::string(" __stdcall")) + " GetModuleData(void) { return &" + std::string(pszSymbolName) + "_DATA; }\n";
	for (size_t idx = 0; idx < pDLL->Exports.size(); idx++)
	{
		if (pDLL->Exports[idx].FunctionNames.size())
			header += "\tconst void* " + pDLL->Exports[idx].FunctionNames.front() + ";\n";
		else
			header += "\tconst void* _" + OutputNumber(idx) + ";\n";
	}
	header += "};\n#pragma pack(pop)\n";

	if (opt_Asm)
	{
		if (!pDLL->is64bit)
		{
			decl += ".MODEL FLAT\n";
			decl += "\n";
			data += "PUBLIC C " + std::string(pszSymbolName) + "_DATA\n";
		}
		else
		{
			data += "PUBLIC " + std::string(pszSymbolName) + "_DATA\n";
		}
		data += "\n";
		data += "_TEXT SEGMENT 'CODE'\n";
		data += "\n";
		data += std::string(pszSymbolName) + "_code PROC C\n";
		data += "\n";
	}
	else
	{
		decl += "#pragma pack(push,1)\nstruct _" + std::string(pszSymbolName) + "_code {\n";
		data += "#pragma section(\".text\")\n";
		data += "extern \"C\" const __declspec(allocate(\".text\")) _" + std::string(pszSymbolName) + "_code " + std::string(pszSymbolName) + "_code =\n{\n";
	}

	// WORD AllocationPages;
	if (opt_Asm)
	{
		data += "\t; WORD AllocationPages;\n";
		data += "\tDW " + OutputWORD((WORD)(pDLL->dwAllocationSize / PAGE_SIZE)) + "\n";
		data += "\n";
	}
	else
	{
		decl += "\tunsigned __int16 " + OutputMemberName() + ";\n";
		data += "\t" + OutputWORD((WORD)(pDLL->dwAllocationSize / PAGE_SIZE)) + ",\n";
	}

	// BYTE SectionsCount;
	if (opt_Asm)
	{
		data += "\t; BYTE SectionsCount;\n";
		data += "\tDB " + OutputBYTE((BYTE)pDLL->Sections.size()) + "\n";
		data += "\n";
	}
	else
	{
		decl += "\tunsigned __int8 " + OutputMemberName() + ";\n";
		data += "\t" + OutputBYTE((BYTE)pDLL->Sections.size()) + ",\n";
	}

	// struct {
	//     WORD SectionBasePage;
	//     BYTE Protection;
	//     WORD SectionPages;
	// } Sections[SectionsCount];
	if (opt_Asm)
	{
		data += "\t; struct {\n";
		data += "\t;     WORD SectionBasePage;\n";
		data += "\t;     BYTE Protection;\n";
		data += "\t;     WORD SectionPages;\n";
		data += "\t; } Sections[" + OutputNumber(pDLL->Sections.size()) + "];\n";
	}
	else
	{
		decl += "\tstruct {\n";
		decl += "\t\tunsigned __int16 " + OutputMemberName() + ";\n";
		decl += "\t\tunsigned __int8 " + OutputMemberName() + ";\n";
		decl += "\t\tunsigned __int16 " + OutputMemberName() + ";\n";
		decl += "\t} " + OutputMemberName() + "[" + OutputNumber(pDLL->Sections.size()) + "];\n";
		data += "\t{\n";
	}
	for (size_t idx = 0; idx < pDLL->Sections.size(); idx++)
	{
		SECTION_TABLE_ENTRY &section = pDLL->Sections[idx];
		if (opt_Asm)
		{
			data += "\tDW " + OutputWORD((WORD)(section.dwSectionBase / PAGE_SIZE)) + "\n";
			data += "\tDB " + OutputBYTE((BYTE)section.dwProtection) + "\n";
			data += "\tDW " + OutputWORD((WORD)(section.dwSectionSize / PAGE_SIZE)) + "\n";
			data += "\n";
		}
		else
		{
			data += "\t\t{ " + OutputWORD((WORD)(section.dwSectionBase / PAGE_SIZE)) + 
					", " + OutputBYTE((BYTE)section.dwProtection) +
					", " + OutputWORD((WORD)(section.dwSectionSize / PAGE_SIZE)) + "}";
			if (idx == pDLL->Sections.size()-1)
				data += "\n";
			else
				data += ",\n";
		}
	}
	if (!opt_Asm)
	{
		data += "\t},\n";
	}

	// BYTE SectionData[SectionDataSize];
	for (std::vector<SECTION_TABLE_ENTRY>::iterator it = pDLL->Sections.begin(); it != pDLL->Sections.end(); it++)
	{
		SECTION_TABLE_ENTRY &section = *it;
		if (opt_Asm)
		{
			data += "\t; BYTE SectionData[" + OutputNumber(section.Data.size()) + "];\n";
			data += OutputData(1, (BYTE*)section.Data.c_str(), section.Data.size());
			data += "\n";
		}
		else
		{
			decl += "\tunsigned __int8 " + OutputMemberName() + "[" + OutputNumber(section.Data.size()) + "];\n";
			data += "\t{\n";
			data += OutputData(2, (BYTE*)section.Data.c_str(), section.Data.size());
			data += "\t},\n";
		}
	}

	// BYTE Relocations[RelocationsSize];
	if (opt_Asm)
	{
		data += "\t; BYTE Relocations[" + OutputNumber(pDLL->Relocations.size()) + "];\n";
		data += OutputData(1, (BYTE*)pDLL->Relocations.c_str(), pDLL->Relocations.size());
		data += "\n";
	}
	else
	{
		decl += "\tunsigned __int8 " + OutputMemberName() + "[" + OutputNumber(pDLL->Relocations.size()) + "];\n";
		data += "\t{\n";
		data += OutputData(2, (BYTE*)pDLL->Relocations.c_str(), pDLL->Relocations.size());
		data += "\t},\n";
	}

	// BYTE  ExportsCount;
	// DWORD Exports[ExportsCount];
	for (size_t idx = 0; idx < pDLL->Exports.size(); idx++)
	{
		if ((idx % 255) == 0)
		{
			if (opt_Asm)
			{
				data += "\t; BYTE ExportsCount;\n";
				data += "\tDB " + OutputBYTE((BYTE)pDLL->Exports.size()) + "\n";
				data += "\n";
				data += "\t; DWORD Exports[" + OutputNumber(pDLL->Exports.size()) + "];\n";
			}
			else
			{
				decl += "\tunsigned __int8 " + OutputMemberName() + ";\n";
				data += "\t" + OutputBYTE((BYTE)pDLL->Exports.size()) + ",\n";
				decl += "\tunsigned __int32 " + OutputMemberName() + "[" + OutputNumber(min(255, pDLL->Exports.size())) + "];\n";
				data += "\t{\n";
			}
		}
		EXPORT_TABLE_ENTRY &exp = pDLL->Exports[idx];
		if (opt_Asm)
		{
			data += "\tDD " + OutputDWORD(exp.dwAddressOfFunction) + " ; ";
		}
		else
		{
			data += "\t\t" + OutputDWORD(exp.dwAddressOfFunction);
			if ((idx % 255) == 255 || idx == pDLL->Exports.size() - 1)
				data += "  //";
			else
				data += ", //";
		}
		if (exp.FunctionNames.size())
		{
			for (std::list<std::string>::iterator it = exp.FunctionNames.begin(); it != exp.FunctionNames.end(); it++)
			{
				if (it != exp.FunctionNames.begin())
					data += ", ";
				data += *it;
			}
			data += "\n";
		}
		else
		{
			data += "@" + OutputNumber(exp.dwExportOrdinal) + "\n";
		}
		if (idx == pDLL->Exports.size() - 1)
		{
			if (opt_Asm)
				data += "\n";
			else
				data += "\t},\n";
		}
	}
	if (opt_Asm)
	{
		data += "\t; BYTE EndOfExports;\n";
		data += "\tDB 0\n";
		data += "\n";
		if (pDLL->is64bit)
		{
			data += "\t; QWORD AddressOfEntryPoint;\n";
			data += "\tDQ " + OutputQWORD(pDLL->pNtHeaders64->OptionalHeader.AddressOfEntryPoint) + "\n";
		}
		else
		{
			data += "\t; DWORD AddressOfEntryPoint;\n";
			data += "\tDD " + OutputDWORD(pDLL->pNtHeaders32->OptionalHeader.AddressOfEntryPoint) + "\n";
		}
		data += "\n";
		data += std::string(pszSymbolName) + "_code ENDP\n";
		data += "\n";
	}
	else
	{
		decl += "\tunsigned __int8 " + OutputMemberName() + ";\n";
		data += "\t0x00,\n";

		if (pDLL->is64bit)
		{
			decl += "\tunsigned __int64 " + OutputMemberName() + ";\n";
			data += "\t" + OutputQWORD(pDLL->pNtHeaders64->OptionalHeader.AddressOfEntryPoint) + "\n";
		}
		else
		{
			decl += "\tunsigned __int32 " + OutputMemberName() + ";\n";
			data += "\t" + OutputDWORD(pDLL->pNtHeaders32->OptionalHeader.AddressOfEntryPoint) + "\n";
		}

		decl += "};\n#pragma pack(pop)\n\n";
		data += "};\n";
	}

	unsigned int lcg = (unsigned int)GetTickCount() ^ (unsigned int)time(NULL);
	unsigned int lcg2 = (unsigned int)GetTickCount() ^ (unsigned int)time(NULL);
	LCG_69069_1(lcg);
	LCG_22695477_1(lcg2);
	lcg = ((LCG_69069_1(lcg) >> 24) & 0xff) |
		  (((LCG_69069_1(lcg) >> 16) & 0xff) << 8) |
		  (((LCG_69069_1(lcg) >> 8) & 0xff) << 16) |
		  ((LCG_69069_1(lcg) & 0xff) << 24);
	if (opt_Asm)
	{
		data += "\n";
		data += std::string(pszSymbolName) + "_DATA PROC C\n";
		data += "\n";
	}
	else
	{
		data += "\nextern \"C\" void _declspec(naked) " + std::string(pszSymbolName) + "_DATA(void)\n{\n\t__asm {\n";
	}
	// random garbage
	if (opt_Asm)
	{
		data += "\tDB " + 
			OutputBYTE((LCG_22695477_1(lcg2) >> 24) & 0xfe) + ", " + // avoid 0xe9
			OutputBYTE((LCG_22695477_1(lcg2) >> 24) & 0xff) + ", " +
			OutputBYTE((LCG_22695477_1(lcg2) >> 24) & 0xff) + ", " +
			OutputBYTE((LCG_22695477_1(lcg2) >> 24) & 0xff) + ", " +
			OutputBYTE((LCG_22695477_1(lcg2) >> 24) & 0xff) + ", " +
			OutputBYTE((LCG_22695477_1(lcg2) >> 24) & 0xff) + ", " +
			OutputBYTE((LCG_22695477_1(lcg2) >> 24) & 0xff) + ", " +
			OutputBYTE((LCG_22695477_1(lcg2) >> 24) & 0xff) + "\n";
	}
	else
	{
		data += "\t\t__emit " + OutputBYTE((LCG_22695477_1(lcg2) >> 24) & 0xfe) + "\n"; // avoid 0xe9
		data += "\t\t__emit " + OutputBYTE((LCG_22695477_1(lcg2) >> 24) & 0xff) + "\n";
		data += "\t\t__emit " + OutputBYTE((LCG_22695477_1(lcg2) >> 24) & 0xff) + "\n";
		data += "\t\t__emit " + OutputBYTE((LCG_22695477_1(lcg2) >> 24) & 0xff) + "\n";
		data += "\t\t__emit " + OutputBYTE((LCG_22695477_1(lcg2) >> 24) & 0xff) + "\n";
		data += "\t\t__emit " + OutputBYTE((LCG_22695477_1(lcg2) >> 24) & 0xff) + "\n";
		data += "\t\t__emit " + OutputBYTE((LCG_22695477_1(lcg2) >> 24) & 0xff) + "\n";
		data += "\t\t__emit " + OutputBYTE((LCG_22695477_1(lcg2) >> 24) & 0xff) + "\n";
	}
	// random seed
	while ((lcg & 0xff) == 0xe9)
		LCG_69069_1(lcg);
	if (opt_Asm)
	{
		data += "\tDB " + 
			OutputBYTE(lcg & 0xff) + ", " +
			OutputBYTE((lcg >> 8) & 0xff) + ", " +
			OutputBYTE((lcg >> 16) & 0xff) + ", " +
			OutputBYTE((lcg >> 24) & 0xff) + "\n";
	}
	else
	{
		data += "\t\t__emit " + OutputBYTE(lcg & 0xff) + "\n";
		data += "\t\t__emit " + OutputBYTE((lcg >> 8) & 0xff) + "\n";
		data += "\t\t__emit " + OutputBYTE((lcg >> 16) & 0xff) + "\n";
		data += "\t\t__emit " + OutputBYTE((lcg >> 24) & 0xff) + "\n";
	}
	// MODULE_code
	LCG_69069_1(lcg);
	for (int c = (int)((lcg >> 24) & 3); c >= 0; c--)
		data += (opt_Asm ? "\tDB " : "\t\t__emit ") + OutputBYTE((BYTE)(LCG_22695477_1(lcg2) >> 16)) + "\n";
	if (opt_Asm)
	{
		data += "\tDB " + OutputBYTE(opcodes[(LCG_22695477_1(lcg2) >> 16) % _countof(opcodes)].v) + "\n";
		data += std::string(opt_64bit ? "\tDQ " : "\tDD ") + std::string(pszSymbolName) + "_code + " + SYMBOL_OFFSET + "\n";
	}
	else
	{
		data += "\t\t" + std::string(opcodes[(LCG_22695477_1(lcg2) >> 16) % _countof(opcodes)].i) + ", far offset " + std::string(pszSymbolName) + "_code + " + SYMBOL_OFFSET + ";\n";
	}
	// VirtualAlloc
	LCG_69069_1(lcg);
	for (int c = (int)((lcg >> 24) & 3); c >= 0; c--)
		data += (opt_Asm ? "\tDB " : "\t\t__emit ") + OutputBYTE((BYTE)(LCG_22695477_1(lcg2) >> 16)) + "\n";
	if (opt_Asm)
	{
		if (opt_64bit)
		{
			decl += "EXTRN __imp_VirtualAlloc:PROC\n";
			data += "\tDB " + OutputBYTE(opcodes[(LCG_22695477_1(lcg2) >> 16) % _countof(opcodes)].v) + "\n";
			data += std::string("\tDQ ") + "__imp_VirtualAlloc + " + SYMBOL_OFFSET + " ; VirtualAlloc [kernel32.dll]\n";
		}
		else
		{
			decl += "EXTRN __imp__VirtualAlloc@16:PROC\n";
			data += "\tDB " + OutputBYTE(opcodes[(LCG_22695477_1(lcg2) >> 16) % _countof(opcodes)].v) + "\n";
			data += std::string("\tDD ") + "__imp__VirtualAlloc@16 + " + SYMBOL_OFFSET + " ; VirtualAlloc [kernel32.dll]\n";
		}
	}
	else
	{
		data += "\t\t" + std::string(opcodes[(LCG_22695477_1(lcg2) >> 16) % _countof(opcodes)].i) + ", far offset VirtualAlloc + " + SYMBOL_OFFSET + ";\n";
	}
	// VirtualProtect
	LCG_69069_1(lcg);
	for (int c = (int)((lcg >> 24) & 3); c >= 0; c--)
		data += (opt_Asm ? "\tDB " : "\t\t__emit ") + OutputBYTE((BYTE)(LCG_22695477_1(lcg2) >> 16)) + "\n";
	if (opt_Asm)
	{
		if (opt_64bit)
		{
			decl += "EXTRN __imp_VirtualProtect:PROC\n";
			data += "\tDB " + OutputBYTE(opcodes[(LCG_22695477_1(lcg2) >> 16) % _countof(opcodes)].v) + "\n";
			data += std::string("\tDQ ") + "__imp_VirtualProtect + " + SYMBOL_OFFSET + " ; VirtualProtect [kernel32.dll]\n";
		}
		else
		{
			decl += "EXTRN __imp__VirtualProtect@16:PROC\n";
			data += "\tDB " + OutputBYTE(opcodes[(LCG_22695477_1(lcg2) >> 16) % _countof(opcodes)].v) + "\n";
			data += std::string("\tDD ") + "__imp__VirtualProtect@16 + " + SYMBOL_OFFSET + " ; VirtualProtect [kernel32.dll]\n";
		}
	}
	else
	{
		data += "\t\t" + std::string(opcodes[(LCG_22695477_1(lcg2) >> 16) % _countof(opcodes)].i) + ", far offset VirtualProtect + " + SYMBOL_OFFSET + ";\n";
	}
	if (opt_Asm)
		data += "\n";

	// dllimports table
	if (opt_Asm)
		data += "\t; DLL imports\n";
	for (IMPORT_BLOCK & import : pDLL->DllImports)
	{
		DWORD dwImportReference = import.dwImportReference;
		// import reference base
		LCG_69069_1(lcg);
		for (int c = (int)((lcg >> 24) & 3); c >= 0; c--)
			data += (opt_Asm ? "\tDB " : "\t\t__emit ") + OutputBYTE((BYTE)(LCG_22695477_1(lcg2) >> 16)) + "\n";
		if (opt_Asm)
		{
			data += "\tDB 03DH\n";
			data += std::string(opt_64bit ? "\tDQ " : "\tDD ") + OutputOFFSET(dwImportReference) + " + " + SYMBOL_OFFSET + "\n";
		}
		else
		{
			data += "\t\tcmp eax, " + OutputOFFSET(dwImportReference) + " + " + SYMBOL_OFFSET + "\n";
		}
		// imported functions
		for (auto it = import.Imports.begin(); it != import.Imports.end(); it++, dwImportReference += sizeof(DWORD))
		{
			IMPORT_TABLE_ENTRY & fnc = *it;
			LCG_69069_1(lcg);
			for (int c = (int)((lcg >> 24) & 3); c >= 0; c--)
				data += (opt_Asm ? "\tDB " : "\t\t__emit ") + OutputBYTE((BYTE)(LCG_22695477_1(lcg2) >> 16)) + "\n";
			if (opt_Asm)
			{
				decl += "EXTRN " + fnc.SymbolName + ":PROC\n";
				data += "\tDB " + OutputBYTE(opcodes[(LCG_22695477_1(lcg2) >> 16) % _countof(opcodes)].v) + "\n";
				data += std::string(opt_64bit ? "\tDQ " : "\tDD ") + fnc.SymbolName + " + " + SYMBOL_OFFSET +
					" ; " + fnc.ImportName + " [" + fnc.ModuleName + "]\n";
			}
			else
			{
				data += "\t\t" + std::string(opcodes[(LCG_22695477_1(lcg2) >> 16) % _countof(opcodes)].i) + ", far offset " + fnc.ImportName + 
//					" + " + SYMBOL_OFFSET +
					"; // " + fnc.ImportName + " [" + fnc.ModuleName + "]\n";
			}
		}
	}
	// end of dllimports
	LCG_69069_1(lcg);
	for (int c = (int)((lcg >> 24) & 3); c >= 0; c--)
		data += (opt_Asm ? "\tDB " : "\t\t__emit ") + OutputBYTE((BYTE)(LCG_22695477_1(lcg2) >> 16)) + "\n";
	if (opt_Asm)
	{
		data += "\tDB 03Dh\n";
		data += std::string(opt_64bit ? "\tDQ " : "\tDD ") + OutputOFFSET(0) + " + " + SYMBOL_OFFSET + "\n";
		data += "\n";
	}
	else
	{
		data += "\t\tcmp eax, 0 + " + SYMBOL_OFFSET + "\n";
	}

	// imports table
	if (opt_Asm)
		data += "\t; imports\n";
	for (IMPORT_BLOCK & import : pDLL->Imports)
	{
		DWORD dwImportReference = import.dwImportReference;
		// import reference base
		LCG_69069_1(lcg);
		for (int c = (int)((lcg >> 24) & 3); c >= 0; c--)
			data += (opt_Asm ? "\tDB " : "\t\t__emit ") + OutputBYTE((BYTE)(LCG_22695477_1(lcg2) >> 16)) + "\n";
		if (opt_Asm)
		{
			data += "\tDB 03DH\n";
			data += std::string(opt_64bit ? "\tDQ " : "\tDD ") + OutputOFFSET(dwImportReference) + " + " + SYMBOL_OFFSET + "\n";
		}
		else
		{
			data += "\t\tcmp eax, " + OutputDWORD(dwImportReference) + " + " + SYMBOL_OFFSET + "\n";
		}
		// imported functions
		for (auto it = import.Imports.begin(); it != import.Imports.end(); it++, dwImportReference += sizeof(DWORD))
		{
			IMPORT_TABLE_ENTRY & fnc = *it;
			LCG_69069_1(lcg);
			for (int c = (int)((lcg >> 24) & 3); c >= 0; c--)
				data += (opt_Asm ? "\tDB " : "\t\t__emit ") + OutputBYTE((BYTE)(LCG_22695477_1(lcg2) >> 16)) + "\n";
			if (opt_Asm)
			{
				decl += "EXTRN " + fnc.SymbolName + ":PROC\n";
				data += "\tDB " + OutputBYTE(opcodes[(LCG_22695477_1(lcg2) >> 16) % _countof(opcodes)].v) + "\n";
				data += std::string(opt_64bit ? "\tDQ " : "\tDD ") + fnc.SymbolName + " + " + SYMBOL_OFFSET +
					" ; " + fnc.ImportName + "\n";
			}
			else
			{
				data += "\t\t" + std::string(opcodes[(LCG_22695477_1(lcg2) >> 16) % _countof(opcodes)].i) + ", far offset " + fnc.ImportName +
//					" + " + SYMBOL_OFFSET +
					"; // " + fnc.ImportName + "\n";
			}
		}
	}
	// end of imports
	LCG_69069_1(lcg);
	for (int c = (int)((lcg >> 24) & 3); c >= 0; c--)
		data += (opt_Asm ? "\tDB " : "\t\t__emit ") + OutputBYTE((BYTE)(LCG_22695477_1(lcg2) >> 16)) + "\n";
	if (opt_Asm)
	{
		data += "\tDB 03DH\n";
		data += std::string(opt_64bit ? "\tDQ " : "\tDD ") + OutputOFFSET(0) + " + " + SYMBOL_OFFSET + "\n";
		data += "\n";
	}
	else
	{
		data += "\t\tcmp eax, 0 + " + SYMBOL_OFFSET + "\n";
	}

	if (opt_Asm)
	{
		decl += "\n";
		data += std::string(pszSymbolName) + "_DATA ENDP\n";
		data += "\n";
		data += "_TEXT ENDS\n";
		data += "\n";
		data += "END\n";
	}
	else
	{
		data += "\t}\n}\n";
	}

	if (_sopen_s(&fd, pszModuleHeaderPath, _O_CREAT | _O_TRUNC | _O_BINARY | _O_WRONLY, _SH_DENYNO, _S_IREAD | _S_IWRITE))
	{
		fprintf(stderr, "error: could not create file \"%s\"\n", pszModuleHeaderPath);
		return false;
	}
	if (_write(fd, header.c_str(), header.size()) != header.size())
	{
		_close(fd);
		fprintf(stderr, "error: could not write to file \"%s\"\n", pszModuleHeaderPath);
		return false;
	}
	_close(fd);

	if (_sopen_s(&fd, pszModuleDataPath, _O_CREAT | _O_TRUNC | _O_BINARY | _O_WRONLY, _SH_DENYNO, _S_IREAD | _S_IWRITE))
	{
		fprintf(stderr, "error: could not create file \"%s\"\n", pszModuleDataPath);
		return false;
	}
	if (_write(fd, decl.c_str(), decl.size()) != decl.size())
	{
		_close(fd);
		fprintf(stderr, "error: could not write to file \"%s\"\n", pszModuleDataPath);
		return false;
	}
	if (_write(fd, data.c_str(), data.size()) != data.size())
	{
		_close(fd);
		fprintf(stderr, "error: could not write to file \"%s\"\n", pszModuleDataPath);
		return false;
	}
	_close(fd);

	return ret;
}

bool ExportSingleSection(DLL *pDLL, const char* pszModuleDataPath, const char* pszSymbolName)
{
	PIMAGE_SECTION_HEADER pSection;
	std::string data;
	int fd;

	if (pDLL->is64bit)
		pSection = (PIMAGE_SECTION_HEADER)IMAGE_FIRST_SECTION(pDLL->pNtHeaders64);
	else
		pSection = (PIMAGE_SECTION_HEADER)IMAGE_FIRST_SECTION(pDLL->pNtHeaders32);

	if (opt_DiBits || opt_LZMA)
		data += "#define " + std::string(pszSymbolName) + "_uncompressed_size " + OutputNumber(pDLL->Sections[0].dwSectionSize) + "\n";
	data += "const unsigned char " + std::string(pszSymbolName) + "[" + OutputNumber(pDLL->Sections[0].Data.size()) + "] =\n{\n";
	data += OutputData(1, (BYTE*)pDLL->Sections[0].Data.data(), pDLL->Sections[0].Data.size());
	data += "};\n";

	if (_sopen_s(&fd, pszModuleDataPath, _O_CREAT | _O_TRUNC | _O_BINARY | _O_WRONLY, _SH_DENYNO, _S_IREAD | _S_IWRITE))
	{
		fprintf(stderr, "error: could not create file \"%s\"\n", pszModuleDataPath);
		return false;
	}
	if (_write(fd, data.c_str(), data.size()) != data.size())
	{
		_close(fd);
		fprintf(stderr, "error: could not write to file \"%s\"\n", pszModuleDataPath);
		return false;
	}
	_close(fd);

	return true;
}


void usage(void)
{
	fprintf(stdout, 
		"usage: dll2mem [options] file.dll [output-header-file (if not single section)] output-data-file symbolname\n"
		"options:\n"
		"   -s | --single-section   - for DLLs with single section, no relocation, no imports and no exports\n"
		"   -d | --dibits           - use di-bits compression (default)\n"
		"   -z | --lzma             - use LZMA compression\n"
		"   -n | --no-compression   - do not compress\n"
		"   -c | --cpp              - CPP output format (default for 32bit)\n"
		"   -a | --asm              - ASM output format (default for 64bit)\n"
		"   -i | --imports obj-file - COFF object file containing import symbols translation table used for ASM output\n"
		"   -32                     - 32-bit output (default)\n"
		"   -64                     - 64-bit output\n"
		"\n"
		);
}


int main(int argc, char* argv[])
{
	int argIdx = 1;
	DLL *pDLL;
	COFF *pCOFF;
	const char *pDllPath = NULL;
	const char *pModuleDataPath = NULL;
	const char *pModuleHeaderPath = NULL;
	const char *pImportsTablePath = NULL;
	const char *pSymbolName = NULL;
	IMPORT_TRANSLATION_TABLE xlat;

	while (argc > argIdx) 
	{ 
		if (!strcmp(argv[argIdx], "-s") || !strcmp(argv[argIdx], "--single-section"))
		{
			if (!opt_SingleSection && pModuleHeaderPath)
				pModuleDataPath = pModuleHeaderPath;
			opt_SingleSection = true;
			++argIdx;
		}
		else if (!strcmp(argv[argIdx], "-z") || !strcmp(argv[argIdx], "--lzma"))
		{
			opt_DiBits = false;
			opt_LZMA = true;
			++argIdx;
		}
		else if (!strcmp(argv[argIdx], "-d") || !strcmp(argv[argIdx], "--dibits"))
		{
			opt_DiBits = true;
			opt_LZMA = false;
			++argIdx;
		}
		else if (!strcmp(argv[argIdx], "-n") || !strcmp(argv[argIdx], "--no-compression"))
		{
			opt_DiBits = false;
			opt_LZMA = false;
			++argIdx;
		}
		else if (!strcmp(argv[argIdx], "-i") || !strcmp(argv[argIdx], "--imports"))
		{
			if (++argIdx >= argc)
			{
				usage();
				exit(1);
			}
			pImportsTablePath = argv[argIdx];
			++argIdx;
		}
		else if (!strcmp(argv[argIdx], "-c") || !strcmp(argv[argIdx], "--cpp"))
		{
			if (!opt_64bit)
				opt_Asm = false;
			argIdx++;
		}
		else if (!strcmp(argv[argIdx], "-a") || !strcmp(argv[argIdx], "--asm"))
		{
			opt_Asm = true;
			++argIdx;
		}
		else if (!strcmp(argv[argIdx], "-64"))
		{
			opt_64bit = true;
			opt_Asm = true;
			++argIdx;
		}
		else if (!strcmp(argv[argIdx], "-32"))
		{
			opt_64bit = false;
			++argIdx;
		}
		else
		{
			if (!pDllPath)
				pDllPath = argv[argIdx];
			else if (!opt_SingleSection && !pModuleHeaderPath)
				pModuleHeaderPath = argv[argIdx];
			else if (!pModuleDataPath)
				pModuleDataPath = argv[argIdx];
			else if (!pSymbolName)
				pSymbolName = argv[argIdx];
			else
			{
				usage();
				exit(1);
			}
			++argIdx;
		}
	}

	if (!pDllPath || (!opt_SingleSection && !pModuleHeaderPath) || !pModuleDataPath || !pSymbolName)
	{
		usage();
		exit(1);
	}

	if (pImportsTablePath)
	{
		pCOFF = LoadCOFFToMemory(pImportsTablePath);
		if (!pCOFF)
			exit(1);
		if (!ParseCOFF(pCOFF))
			exit(1);
		if (!BuildImportsTranslationTable(pCOFF, &xlat))
			exit(1);
	}

	pDLL = LoadDLLToMemory(pDllPath);
	if (!pDLL)
		exit(1);
	if (!ParseDll(pDLL, pImportsTablePath ? &xlat : nullptr))
		exit(1);

	if (opt_SingleSection)
	{
		if (!ExportSingleSection(pDLL, pModuleDataPath, pSymbolName))
			exit(1);
	}
	else
	{
		if (!ExportDLL(pDLL, pModuleHeaderPath, pModuleDataPath, pSymbolName))
			exit(1);
	}

	return 0;
}
