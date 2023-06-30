#include "file.h"
#include "data_load.h"
#include "data_storage.h"

bool manageFile(dataType  ** imageDataPtr, const size_t length, const size_t width,
	const size_t height, unsigned char * pathPtr, OperationType operation, LoadDataType dType, Storage_Flags flags)
{
	bool status = false; // Initial Status, only changed to true if the operation is successful
    switch (operation)
    {
        case LOAD_DATA_RAW:
            status = load3dDataArrayRAW(imageDataPtr, length, width, height, pathPtr, dType);
            break;
        case LOAD_2D_DATA_PGM:
            switch (dType)
            {
                case BINARY_DATA:
                    //TODO:I want to be able to read 2D image (PGM raw, RAW)
                    break;
                case ASCII_DATA:
                    status = load2dPGMASCII(imageDataPtr, width, height, pathPtr);
                    break;
                default:
                    break;
            }
            break;
        case STORE_DATA_RAW:
            switch (dType)
            {
                case BINARY_DATA:
                    status = store3dDataArrayD(imageDataPtr, length, width, height, pathPtr, flags);
                    break;
                case ASCII_DATA:
                    status = store3dDataArrayASCII(imageDataPtr, length, width, height, pathPtr, flags);
                    break;
                default:
                    break;
            }
            break;
        case STORE_2D_DATA_PGM:
            switch (dType)
            {
            case BINARY_DATA:
                //TODO: I want to be able to store 2D image (PGM raw, RAW)
                break;
            case ASCII_DATA:
                status = store2dPGMASCII(imageDataPtr, width, height, pathPtr);
                break;
            default:
                break;
            }
            break;
        default:
            break;
    }
	return status;
}

void convertTodataType(unsigned char ** dataPtrUC, dataType ** dataPtrD, const size_t dimXY, const size_t height)
{
	for (size_t k = 0; k < height; k++)
	{
		for (size_t i = 0; i < dimXY; i++)
		{
			dataPtrD[k][i] = (dataType)dataPtrUC[k][i];
		}
	}
}

void convertFromShortTodatatype(short** dataPtrSH, dataType** dataPtrF, const size_t dimXY, const size_t height)
{
	for (size_t k = 0; k < height; k++)
	{
		for (size_t i = 0; i < dimXY; i++)
		{
			dataPtrF[k][i] = (dataType)dataPtrSH[k][i];
		}
	}
}
