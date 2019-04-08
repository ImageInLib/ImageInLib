#include "headers.h"
#include "check.h"
#include "rcfile.h"
#include "typeinfo.h"
#include <shlwapi.h>
#include <io.h>


std::wstring unquoteString(const wchar_t *str)
{
	std::wstring s;
	if (*str == '"')
	{
		str++;
		while (*str)
		{
			if (*str == '\\' && *(str+1) == '"')
			{
				s += L"\"";
				str++;
			} 
			else if (*str == '\\' && *(str+1) == '\\')
			{
				s += L"\\";
				str++;
			}
			else if (*str == '"' && *(str+1) == '-')
			{
				s += L"'-";
				str++;
			}
//			else if (*str == '"' && *(str+1) >= '0' && *(str+1) <= '9' )
//			{
//				s += '\''; s+= *(str+1);
//				str++;
//			}
			else if (*str == '"' && *(str+1) == '"')
			{
				s += L"\"";
				str++;
			}
			else if (*str == '"' && *(str+1) == 0)
			{
				;
			}
			else
			{
				s += *str;
			}
			str++;
		}
	}
	else
	{
		s = str;
	}
	return s;
}

static int checkBitmap(const wchar_t *pszBitmapName, const wchar_t *pszFilePath)
{
	BITMAPFILEHEADER bmpFileHeader;
	BITMAPINFOHEADER bmpInfoHeader;
	DWORD bmpBitFields[3];
	int fd;

	//
	// Open source bitmap file
	//
	if (_wsopen_s(&fd, pszFilePath, _O_RDONLY | _O_BINARY, _SH_DENYNO, _S_IREAD) != 0 || fd == -1)
	{
		wprintf(L"could not open bitmap %s in file '%s'\n", pszBitmapName, pszFilePath);
		return 1;
	}

	//
	// Read source bitmap file headers
	//
	if (sizeof(BITMAPFILEHEADER) != _read(fd, &bmpFileHeader, sizeof(bmpFileHeader)))
	{
		wprintf(L"could not read BITMAPFILEHEADER from bitmap %s in file '%s'\n", pszBitmapName, pszFilePath);
		_close(fd);
		return 1;
	}
	if (bmpFileHeader.bfType != 'MB')
	{
		wprintf(L"unknown format of bitmap %s in file '%s'\n", pszBitmapName, pszFilePath);
		_close(fd);
		return 1;
	}
	if (sizeof(bmpInfoHeader.biSize) != _read(fd, &bmpInfoHeader.biSize, sizeof(bmpInfoHeader.biSize)))
	{
		wprintf(L"could not read BITMAPINFOHEADER from bitmap %s in file '%s'\n", pszBitmapName, pszFilePath);
		_close(fd);
		return 1;
	}
	if (bmpInfoHeader.biSize < sizeof(BITMAPINFOHEADER))
	{
		wprintf(L"unknown format of bitmap %s in file '%s'\n", pszBitmapName, pszFilePath);
		_close(fd);
		return 1;
	}
	if (sizeof(BITMAPINFOHEADER) - sizeof(bmpInfoHeader.biSize) != 
		_read(fd, &bmpInfoHeader.biWidth, sizeof(BITMAPINFOHEADER) - sizeof(bmpInfoHeader.biSize)))
	{
		wprintf(L"could not read BITMAPINFOHEADER from bitmap %s in file '%s'\n", pszBitmapName, pszFilePath);
		_close(fd);
		return 1;
	}
	if (bmpInfoHeader.biBitCount != 32)
	{
		_close(fd);
		return 0;
	}
	if (bmpInfoHeader.biCompression != BI_BITFIELDS)
	{
		wprintf(L"bitmap %s does not have BI_BITFIELDS format\n", pszBitmapName);
		_close(fd);
		return 1;
	}
	if (bmpInfoHeader.biCompression == BI_BITFIELDS)
	{
		if (sizeof(BITMAPFILEHEADER) + bmpInfoHeader.biSize != _lseek(fd, sizeof(BITMAPFILEHEADER) + bmpInfoHeader.biSize, SEEK_SET))
		{
			wprintf(L"could not read color masks from bitmap %s in file '%s'\n", pszBitmapName, pszFilePath);
			_close(fd);
			return 1;
		}
		if (3*sizeof(DWORD) != _read(fd, &bmpBitFields, sizeof(bmpBitFields)))
		{
			wprintf(L"could not read color masks from bitmap %s in file '%s'\n", pszBitmapName, pszFilePath);
			_close(fd);
			return 1;
		}
		if (bmpBitFields[0] != 0x00ff0000 ||
			bmpBitFields[1] != 0x0000ff00 ||
			bmpBitFields[2] != 0x000000ff)
		{
			wprintf(L"bitmap %s uses non standard color masks\n", pszBitmapName);
			_close(fd);
			return 1;
		}
	}

	_close(fd);
	return 0;
}

int checkBitmaps(RCFile_C *pRC, const wchar_t *pszBaseDir)
{
	int res = 0;

	for (int resIdx = 0; resIdx < pRC->GetResourceCount(); resIdx++)
	{
		Resource_C *pRes;
		pRes = pRC->GetResource(resIdx);
		if (dynamic_cast<Bitmap_C*>(pRes))
		{
			Bitmap_C *pBitmap = dynamic_cast<Bitmap_C*>(pRes);
			std::wstring filePath;
			filePath = pszBaseDir;
			if (filePath.size() && filePath.at(filePath.size()-1) != '\\')
				filePath += L"\\";
			filePath += unquoteString(pBitmap->GetFileName());
			if (checkBitmap(pBitmap->GetName(), filePath.c_str()))
				res = 1;
		}
	}

	return res;
}


#pragma pack(push)
#pragma pack(1)

typedef struct CURSORDIR {
	WORD reserved;          // should always be 0
	WORD nImageType;        // 1 for icon (.ICO) image, 2 for cursor (.CUR) image
	WORD nImagesCount;		// number of images in file
} CURSORDIR;

typedef struct CURSORDIRENTRY {
	BYTE nWidth;        // image width in pixels (value 0 means image width is 256 pixels)
	BYTE nHeight;       // image height in pixels (value 0 means image height is 256 pixels)
	BYTE nColors;       // number of colors in the color palette, should be 0 if the image does not use a color palette
	BYTE reserved;      // should always be 0
	WORD nHotSpotX;     // the horizontal coordinates of the hotspot in number of pixels from the left
	WORD nHotSpotY;     // the vertical coordinates of the hotspot in number of pixels from the top
	DWORD nImageSize;   // the size of the image's data in bytes
	DWORD nImageOffset; // the offset of BMP or PNG data from the beginning of the CUR file
} CURSORDIRENTRY;

#pragma pack(pop)

static int checkCursor(const wchar_t *pszCursorName, const wchar_t *pszFilePath)
{
	CURSORDIR cursorDir;
	int fd;

	//
	// Open source bitmap file
	//
	if (_wsopen_s(&fd, pszFilePath, _O_RDONLY | _O_BINARY, _SH_DENYNO, _S_IREAD) != 0 || fd == -1)
	{
		wprintf(L"could not open cursor %s in file '%s'\n", pszCursorName, pszFilePath);
		return 1;
	}

	//
	// Read source bitmap file headers
	//
	if (sizeof(CURSORDIR) != _read(fd, &cursorDir, sizeof(cursorDir)))
	{
		wprintf(L"could not read CURSORDIR from cursor %s in file '%s'\n", pszCursorName, pszFilePath);
		_close(fd);
		return 1;
	}
	if (cursorDir.reserved != 0 || cursorDir.nImageType != 2)
	{
		wprintf(L"unknown format of cursor %s in file '%s'\n", pszCursorName, pszFilePath);
		return 1;
	}
	for (int imgIdx = 0; imgIdx < cursorDir.nImagesCount; imgIdx++)
	{
		CURSORDIRENTRY cursorDirEntry;
		BITMAPINFOHEADER bmpInfoHeader;

		if (_lseek(fd, sizeof(CURSORDIR) + sizeof(CURSORDIRENTRY) * imgIdx, SEEK_SET) != sizeof(CURSORDIR) + sizeof(CURSORDIRENTRY) * imgIdx)
		{
			wprintf(L"could not read CURSORDIRENTRY from cursor %s in file '%s'\n", pszCursorName, pszFilePath);
			_close(fd);
			return 1;
		}
		if (sizeof(CURSORDIRENTRY) != _read(fd, &cursorDirEntry, sizeof(CURSORDIRENTRY)))
		{
			wprintf(L"could not read CURSORDIRENTRY from cursor %s in file '%s'\n", pszCursorName, pszFilePath);
			_close(fd);
			return 1;
		}
		if (cursorDirEntry.reserved != 0)
		{
			wprintf(L"unknown format of cursor %s in file '%s'\n", pszCursorName, pszFilePath);
			_close(fd);
			return 1;
		}

		if (_lseek(fd, cursorDirEntry.nImageOffset, SEEK_SET) != cursorDirEntry.nImageOffset)
		{
			wprintf(L"could not read BITMAPINFOHEADER from cursor %s in file '%s'\n", pszCursorName, pszFilePath);
			_close(fd);
			return 1;
		}

		if (sizeof(BITMAPINFOHEADER) != _read(fd, &bmpInfoHeader, sizeof(BITMAPINFOHEADER)))
		{
			wprintf(L"could not read BITMAPINFOHEADER from bitmap %s in file '%s'\n", pszCursorName, pszFilePath);
			_close(fd);
			return 1;
		}

		if (bmpInfoHeader.biBitCount != 32 && bmpInfoHeader.biBitCount != 1)
		{
			wprintf(L"cursor %s does not have 32-bit or 1-bit format\n", pszCursorName);
			_close(fd);
			return 1;
		}
	}

	_close(fd);
	return 0;
}

int checkCursors(RCFile_C *pRC, const wchar_t *pszBaseDir)
{
	int res = 0;

	for (int resIdx = 0; resIdx < pRC->GetResourceCount(); resIdx++)
	{
		Resource_C *pRes;
		pRes = pRC->GetResource(resIdx);
		if (dynamic_cast<Cursor_C*>(pRes))
		{
			Cursor_C *pCursor = dynamic_cast<Cursor_C*>(pRes);
			std::wstring filePath;
			filePath = pszBaseDir;
			if (filePath.size() && filePath.at(filePath.size()-1) != '\\')
				filePath += L"\\";
			filePath += unquoteString(pCursor->GetFileName());
			if (checkCursor(pCursor->GetName(), filePath.c_str()))
				res = 1;
		}
	}

	return res;
}
