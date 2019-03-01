#include "vtk_params.h"
#include "file.h"
#include "data_load.h"
#include <stdbool.h>
#include "data_storage.h"

bool manageFile(void  ** imageDataPtr, const size_t length, const size_t width,
	const size_t height, unsigned char * pathPtr, VTKHeaderLines *lines, operationType operation, loadDataType dType)
{
	bool status = false; // Initial Status, only changed to true if the operation is successful
	switch (operation)
	{
	case LOAD_DATA_VTK:
		status = load3dDataArrayVTK((unsigned char **)imageDataPtr, length, width, height, pathPtr, lines);
		break;
	case LOAD_DATA_RAW:
		status = load3dDataArrayRAW(imageDataPtr, length, width, height, pathPtr, dType);
		break;
	case STORE_DATA_VTK:
		status = store3dRealDataVtkD((double **)imageDataPtr, length, width, height, pathPtr, lines);
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