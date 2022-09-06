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

void rescaleNewRange(dataType** imageDataPtr, size_t imageLength, size_t imageWidth, size_t imageHeight, dataType minNew, dataType maxNew) {
	size_t k, i, j, xd;
	dataType maxData = -10000; 
	dataType minData = 10000;

	// Find the Max Intensity
	for (k = 0; k < imageHeight; k++) {
		for (i = 0; i < imageLength; i++) {
			for (j = 0; j < imageWidth; j++) {

				// 1D Conversion of row and column
				xd = x_new(i, j, imageLength);
				if (imageDataPtr[k][xd] > maxData) {
					maxData = imageDataPtr[k][xd];
				}
				if (imageDataPtr[k][xd] < minData) {
					minData = imageDataPtr[k][xd];
				}
			}
		}
	}
	// Rescale from min_new to max_new
	dataType diffOld = maxData - minData;
	dataType diffNew = maxNew - minNew;
	dataType scale_factor = (diffNew) / (diffOld);
	for (k = 0; k < imageHeight; k++) {
		for (i = 0; i < imageLength; i++) {
			for (j = 0; j < imageWidth; j++) {
				// 1D Conversion of row and column
				xd = x_new(i, j, imageLength);
				imageDataPtr[k][xd] = scale_factor * (imageDataPtr[k][xd] - maxData) + maxNew;
				// Alternatively
				//dta[k][xd] = scale_factor * (dta[k][xd] - min_dta) + min_new; }
			}
		}
	}
}
