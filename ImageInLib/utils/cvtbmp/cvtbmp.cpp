// cvtbmp.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <string>

BYTE buffer[65536];

enum CONVERSION_TYPE {
	CVT_UNSPECIFIED,
	CVT_TO_BITFIELDS,
	CVT_TO_RGB
};

enum EXIT_CODE {
	EXIT_NORMAL = 0,
	EXIT_INVALID_OPTIONS = 1,	
	EXIT_SOURCE_OPEN_ERROR = 2,
	EXIT_SOURCE_READ_ERROR = 3,
	EXIT_SOURCE_FORMAT_ERROR = 4,
	EXIT_DEST_CREATE_ERROR = 5,
	EXIT_DEST_WRITE_ERROR = 6,
	EXIT_UNSPECIFIED_ERROR = 99
};

void usage(void)
{
	fprintf(stdout,
"Usage: cvtbmp [options] SOURCE DEST\n"
"Changes format of BMP file from BI_RGB to BI_BITFIELDS and vice versa.\n"
"\n"
"Conversion options:\n"
"  --bitfields   converts to BI_BITFIELDS format (default)\n"
"  --rgb         converts to BI_RGB format\n"
//"  --target-directory DIRECTORY\n"
//"                creates DEST in DIRECTORY\n"
"  --help        display this help and exit\n"
"\n"
"Exit codes:\n"
"  0  Normal, no errors or warnings detected.\n"
"  1  Invalid options were specified on the command line.\n"
"  2  Could not open SOURCE file\n"
"  3  Could not read SOURCE file.\n"
"  4  SOURCE file has invalid format.\n"
"  5  Could not create DEST file.\n"
"  6  Could not write to DEST file.\n"
"  99 Unspecified error encountered. \n"
"\n"
"Notes:\n"
"  Only BMP files with 32 bits per pixel are supported.\n"
"\n"
	);
}

int wmain(int argc, WCHAR* argv[])
{
	CONVERSION_TYPE convType = CVT_UNSPECIFIED;
	const WCHAR *sourceFile = NULL;
	const WCHAR *destFile = NULL;
	const WCHAR *targetDirectory = NULL;
	std::wstring destPath;
	int fdSource;
	int fdDest;
	BITMAPFILEHEADER bmpFileHeader;
	BITMAPINFOHEADER bmpInfoHeader;
	DWORD bmpBitFields[3];

	//
	// Parse command line
	//
	while (argc > 1)
	{
		if (wcscmp(argv[1], L"--help") == 0)
		{
			usage();
			exit(EXIT_NORMAL);
		}
		else if (wcscmp(argv[1], L"--bitfields") == 0)
		{
			if (convType != CVT_UNSPECIFIED && convType != CVT_TO_BITFIELDS)
			{
				fprintf(stderr,	"cvtbmp: multiple conversion options specified\n");
				exit(EXIT_INVALID_OPTIONS);
			}
			convType = CVT_TO_BITFIELDS;
		}
		else if (wcscmp(argv[1], L"--rgb") == 0)
		{
			if (convType != CVT_UNSPECIFIED && convType != CVT_TO_RGB)
			{
				fprintf(stderr,	"cvtbmp: multiple conversion options specified\n");
				exit(EXIT_INVALID_OPTIONS);
			}
			convType = CVT_TO_RGB;
		}
		else if (wcscmp(argv[1], L"--target-directory") == 0)
		{
			if (argc < 3)
			{
				fprintf(stderr,	"cvtbmp: missing target directory\n");
				exit(EXIT_INVALID_OPTIONS);
			}
			targetDirectory = argv[2];
			argc--;
			argv++;
		}
		else
		{
			break;
		}
		argc--;
		argv++;
	}
	if (argc != 3)
	{
		fprintf(stderr,	"cvtbmp: missing file arguments\n"
		                "Try `cvtbmp --help' for more information.\n"
			   );
		exit(EXIT_INVALID_OPTIONS);
	}
	sourceFile = argv[1];
	destFile = argv[2];
	if (targetDirectory != NULL && PathIsRelative(destFile))
	{
		PathRemoveBackslash((LPWSTR)targetDirectory);
		destPath = targetDirectory;
		destPath += '\\';
		destPath += PathFindFileName(destFile);
		destFile = destPath.c_str();
	}

	//
	// Open source bitmap file
	//
	if (_wsopen_s(&fdSource, sourceFile, _O_RDONLY | _O_BINARY, _SH_DENYNO, _S_IREAD) != 0 || fdSource == -1)
	{
		fprintf(stderr,	"cvtbmp: could not open source file\n");
		exit(EXIT_SOURCE_OPEN_ERROR);
	}

	//
	// Read source bitmap file headers
	//
	if (sizeof(BITMAPFILEHEADER) != _read(fdSource, &bmpFileHeader, sizeof(bmpFileHeader)))
	{
		fprintf(stderr,	"cvtbmp: could not read BITMAPFILEHEADER from source file\n");
		exit(EXIT_SOURCE_READ_ERROR);
	}
	if (bmpFileHeader.bfType != 'MB')
	{
		fprintf(stderr,	"cvtbmp: invalid source file format, only 'BM' file type is supported\n");
		exit(EXIT_SOURCE_FORMAT_ERROR);
	}
	if (sizeof(bmpInfoHeader.biSize) != _read(fdSource, &bmpInfoHeader.biSize, sizeof(bmpInfoHeader.biSize)))
	{
		fprintf(stderr,	"cvtbmp: could not read BITMAPINFOHEADER from source file\n");
		exit(EXIT_SOURCE_READ_ERROR);
	}
	if (bmpInfoHeader.biSize < sizeof(BITMAPINFOHEADER))
	{
		fprintf(stderr,	"cvtbmp: invalid source file format, file with BITMAPINFOHEADER is required\n");
		exit(EXIT_SOURCE_FORMAT_ERROR);
	}
	if (sizeof(BITMAPINFOHEADER) - sizeof(bmpInfoHeader.biSize) != 
		_read(fdSource, &bmpInfoHeader.biWidth, sizeof(BITMAPINFOHEADER) - sizeof(bmpInfoHeader.biSize)))
	{
		fprintf(stderr,	"cvtbmp: could not read BITMAPINFOHEADER from source file\n");
		exit(EXIT_SOURCE_READ_ERROR);
	}
	if (bmpInfoHeader.biBitCount != 32)
	{
		fprintf(stderr,	"cvtbmp: invalid source file format, only 32bpp file is supported\n");
		exit(EXIT_SOURCE_FORMAT_ERROR);
	}
	if (bmpInfoHeader.biCompression != BI_RGB && bmpInfoHeader.biCompression != BI_BITFIELDS)
	{
		fprintf(stderr,	"cvtbmp: invalid source file format, only BI_RGB and BI_BITFIELDS is supported\n");
		exit(EXIT_SOURCE_FORMAT_ERROR);
	}
	if (bmpInfoHeader.biCompression == BI_BITFIELDS)
	{
		if (sizeof(BITMAPFILEHEADER) + bmpInfoHeader.biSize != _lseek(fdSource, sizeof(BITMAPFILEHEADER) + bmpInfoHeader.biSize, SEEK_SET))
		{
			fprintf(stderr,	"cvtbmp: could not read color masks from source file\n");
			exit(EXIT_SOURCE_READ_ERROR);
		}
		if (3*sizeof(DWORD) != _read(fdSource, &bmpBitFields, sizeof(bmpBitFields)))
		{
			fprintf(stderr,	"cvtbmp: could not read color masks from source file\n");
			exit(EXIT_SOURCE_READ_ERROR);
		}
	}
	if (bmpFileHeader.bfOffBits != _lseek(fdSource, bmpFileHeader.bfOffBits, SEEK_SET))
	{
		fprintf(stderr,	"cvtbmp: could not read bitmap bits from source file\n");
		exit(EXIT_SOURCE_READ_ERROR);
	}

	//
	// Convert file headers
	//
	if (convType == CVT_TO_RGB)
	{
		DWORD dwBitsSize = bmpFileHeader.bfSize - bmpFileHeader.bfOffBits;
		bmpFileHeader.bfOffBits = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER);
		bmpFileHeader.bfSize = bmpFileHeader.bfOffBits + dwBitsSize;
		bmpInfoHeader.biSize = sizeof(BITMAPINFOHEADER);
		bmpInfoHeader.biCompression = BI_RGB;
		bmpInfoHeader.biClrUsed = 0;
		bmpInfoHeader.biClrImportant = 0;
	}
	else
	{
		DWORD dwBitsSize = bmpFileHeader.bfSize - bmpFileHeader.bfOffBits;
		bmpFileHeader.bfOffBits = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER) + sizeof(bmpBitFields);
		bmpFileHeader.bfSize = bmpFileHeader.bfOffBits + dwBitsSize;
		if (bmpInfoHeader.biCompression != BI_BITFIELDS)
		{
			bmpBitFields[0] = 0x00FF0000;
			bmpBitFields[1] = 0x0000FF00;
			bmpBitFields[2] = 0x000000FF;
		}
		bmpInfoHeader.biSize = sizeof(BITMAPINFOHEADER);
		bmpInfoHeader.biCompression = BI_BITFIELDS;
		bmpInfoHeader.biClrUsed = 0;
		bmpInfoHeader.biClrImportant = 0;
	}

	//
	// Create destination file
	//
	if (_wsopen_s(&fdDest, destFile, _O_CREAT | _O_TRUNC | _O_RDWR | _O_BINARY, _SH_DENYNO, _S_IREAD | _S_IWRITE) != 0 || fdSource == -1)
	{
		fprintf(stderr,	"cvtbmp: could not create destination file\n");
		exit(EXIT_DEST_CREATE_ERROR);
	}

	//
	// Write bitmap headers to the destination file
	//
	if (sizeof(BITMAPFILEHEADER) != _write(fdDest, &bmpFileHeader, sizeof(BITMAPFILEHEADER)))
	{
		fprintf(stderr,	"cvtbmp: could not write BITMAPFILEHEADER to destination file\n");
		exit(EXIT_DEST_WRITE_ERROR);
	}
	if (sizeof(BITMAPINFOHEADER) != _write(fdDest, &bmpInfoHeader, sizeof(BITMAPINFOHEADER)))
	{
		fprintf(stderr,	"cvtbmp: could not write BITMAPINFOHEADER to destination file\n");
		exit(EXIT_DEST_WRITE_ERROR);
	}
	if (convType != CVT_TO_RGB && 3*sizeof(DWORD) != _write(fdDest, &bmpBitFields, sizeof(bmpBitFields)))
	{
		fprintf(stderr,	"cvtbmp: could not write color masks to destination file\n");
		exit(EXIT_DEST_WRITE_ERROR);
	}

	//
	// Copy bitmap bits to the destination file
	//
	for (size_t bitmapSize = bmpFileHeader.bfSize - bmpFileHeader.bfOffBits, bytesToRead, bytesRead, bytesWritten;
		 0 != (bytesToRead = min(sizeof(buffer), bitmapSize));
		 bitmapSize -= bytesWritten)
	{
		bytesRead = _read(fdSource, buffer, bytesToRead);
		if (bytesRead != bytesToRead)
		{
			fprintf(stderr,	"cvtbmp: could not read bitmap bits from source file\n");
			exit(EXIT_SOURCE_READ_ERROR);
		}
		bytesWritten = _write(fdDest, buffer, bytesRead);
		if (bytesWritten != bytesRead)
		{
			fprintf(stderr,	"cvtbmp: could not write bitmap bits to destination file\n");
			exit(EXIT_DEST_WRITE_ERROR);
		}
	}

	_close(fdSource);
	_close(fdDest);

	return EXIT_NORMAL;
}

