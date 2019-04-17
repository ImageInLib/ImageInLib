#include "common_vtk.h"
//==============================================================================
// VTK smart pointer
#include "vtkSmartPointer.h"
// VTK Data Reader
#include "vtkGenericDataObjectReader.h"
// VTK Data writer
#include "vtkGenericDataObjectWriter.h"

#include "vtkDataSet.h"
//==============================================================================
void fillPtr(double * ptr, Vtk_File_Info * vtkInfo);
void fillPtr(int * ptr, Vtk_File_Info * vtkInfo);
void fillPtr(float * ptr, Vtk_File_Info * vtkInfo);
void fillPtr(unsigned int * ptr, Vtk_File_Info * vtkInfo);
void fillPtr(unsigned char * ptr, Vtk_File_Info * vtkInfo);
//==============================================================================
int readVtkFile(const char * inputFilePath, Vtk_File_Info * vtkMetaInfo)
{
	// Set up the reader obj
	vtkSmartPointer<vtkGenericDataObjectReader> readGenericVtk = vtkSmartPointer<vtkGenericDataObjectReader>::New();
	// Set the input details
	readGenericVtk->SetFileName(inputFilePath);
	// Read the vtk file
	try
	{
		readGenericVtk->Update();
	}
	catch (const std::exception &)
	{
		std::cout << "Error reading file" << std::endl;
		return EXIT_FAILURE;
	}
	// If successful, save the info to vtkMetaInfo
	vtkDataObject * dataObject = readGenericVtk->GetOutput();
	vtkDataSet * dataSet = vtkDataSet::SafeDownCast(dataObject);
	vtkImageData* imageData = vtkImageData::SafeDownCast(dataSet);
	// Extract the Meta info and store to the vtkMetaInfo
	imageData->GetDimensions(vtkMetaInfo->dimensions);
	imageData->GetSpacing(vtkMetaInfo->spacing);
	imageData->GetOrigin(vtkMetaInfo->origin);
	// Sets the operation
	vtkMetaInfo->operation = copyTo;
	// Assign the data type and pointer to fill
	switch (imageData->GetScalarType())
	{
	case VTK_DOUBLE:
	{
		vtkMetaInfo->vDataType = dta_Dbl;
		// Assign the Data Pointer
		double * ptr_d = (double *)imageData->GetScalarPointer();
		fillPtr(ptr_d, vtkMetaInfo);
		break;
	}
	case VTK_FLOAT:
	{
		vtkMetaInfo->vDataType = dta_Flt;
		// Assign the Data Pointer
		float * ptr_f = (float *)imageData->GetScalarPointer();
		fillPtr(ptr_f, vtkMetaInfo);
		break;
	}
	case VTK_INT:
	{
		vtkMetaInfo->vDataType = dta_Int;
		// Assign the Data Pointer
		int * ptr_i = (int *)imageData->GetScalarPointer();
		fillPtr(ptr_i, vtkMetaInfo);
		break;
	}
	case VTK_UNSIGNED_CHAR:
	{
		vtkMetaInfo->vDataType = dta_UChar;
		// Assign the Data Pointer
		unsigned char * ptr_uc = (unsigned char *)imageData->GetScalarPointer();
		fillPtr(ptr_uc, vtkMetaInfo);
		break;
	}
	case VTK_UNSIGNED_INT:
	{
		vtkMetaInfo->vDataType = dta_UInt;
		// Assign the Data Pointer
		unsigned int * ptr_ui = (unsigned int *)imageData->GetScalarPointer();
		fillPtr(ptr_ui, vtkMetaInfo);
		break;
	}
	default:
		break;
	}
	return EXIT_SUCCESS;
}
int storeVtkFile(const char * outputFilePath, Vtk_File_Info * vtkMetaInfo)
{
	// Creates A new imageData Pointer
	vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
	// FIlls the correct meta data
	createVtkImageData(imageData, vtkMetaInfo);
	// VTK Reader
	// Create a Writer
	vtkSmartPointer<vtkGenericDataObjectWriter> writeGenericVtk = vtkSmartPointer<vtkGenericDataObjectWriter>::New();

	writeGenericVtk->SetInputData(imageData);
	writeGenericVtk->SetFileName(outputFilePath);
	writeGenericVtk->SetFileTypeToASCII();
	try
	{
		writeGenericVtk->Update();
	}
	catch (const std::exception&)
	{
		std::cout << "Error writing file" << std::endl;
		return EXIT_FAILURE;
	}
	return 0;
}
//==============================================================================
int createVtkImageData(vtkImageData * imageData, Vtk_File_Info * vtkInfo)
{
	imageData->SetDimensions(vtkInfo->dimensions);
	imageData->SetSpacing(vtkInfo->spacing);
	imageData->SetOrigin(vtkInfo->origin);
	// Sets the operation
	vtkInfo->operation = copyFrom;
	// Define Type of VTK data type method to use
	switch (vtkInfo->vDataType)
	{
	case dta_Dbl:
	{
		imageData->AllocateScalars(VTK_DOUBLE, 1);
		// Assign the Data Pointer
		double * ptr = (double *)imageData->GetScalarPointer();
		// Fill in the data
		fillPtr(ptr, vtkInfo);
		break;
	}
	case dta_Flt:
	{
		imageData->AllocateScalars(VTK_FLOAT, 1);
		// Assign the Data Pointer
		float * ptr = (float *)imageData->GetScalarPointer();
		// Fill in the data
		fillPtr(ptr, vtkInfo);
		break;
	}
	case dta_Int:
	{
		imageData->AllocateScalars(VTK_INT, 1);
		// Assign the Data Pointer
		int * ptr = (int *)imageData->GetScalarPointer();
		// Fill in the data
		fillPtr(ptr, vtkInfo);
		break;
	}
	case dta_UChar:
	{
		imageData->AllocateScalars(VTK_UNSIGNED_CHAR, 1);
		// Assign the Data Pointer
		unsigned char * ptr = (unsigned char *)imageData->GetScalarPointer();
		// Fill in the data
		fillPtr(ptr, vtkInfo);
		break;
	}
	case dta_UInt:
	{
		imageData->AllocateScalars(VTK_UNSIGNED_INT, 1);
		// Assign the Data Pointer
		unsigned int * ptr = (unsigned int *)imageData->GetScalarPointer();
		// Fill in the data
		fillPtr(ptr, vtkInfo);
		break;
	}
	default:
		break;
	}
	return EXIT_SUCCESS;
}
//==============================================================================
void fillPtr(double * ptr, Vtk_File_Info * vtkInfo)
{
	size_t k, i, j, xd, xyd;
	int height = vtkInfo->dimensions[2], length = vtkInfo->dimensions[1], width = vtkInfo->dimensions[0];
	if (vtkInfo->operation == copyTo)
	{
		// Initialize the pointer to which it is to be copied to
		vtkInfo->dataPointer = (dataType **)malloc(sizeof(dataType *) * height);
		for (i = 0; i < height; i++)
		{
			vtkInfo->dataPointer[i] = (dataType *)malloc(sizeof(dataType) * (length * width));
		}
	}
	for (k = 0; k < height; k++)
	{
		for (i = 0; i < length; i++)
		{
			for (j = 0; j < width; j++)
			{
				xd = x_new(i, j, length);
				xyd = x_flat(i, j, k, length, width);
				if (vtkInfo->operation == copyFrom)
				{
					ptr[xyd] = static_cast<double>(vtkInfo->dataPointer[k][xd]);
				}
				else if (vtkInfo->operation == copyTo)
				{
					vtkInfo->dataPointer[k][xd] = static_cast<dataType>(ptr[xyd]);
				}
			}
		}
	}
}
//==============================================================================
void fillPtr(int * ptr, Vtk_File_Info * vtkInfo)
{
	size_t k, i, j, xd, xyd;
	int height = vtkInfo->dimensions[2], length = vtkInfo->dimensions[1], width = vtkInfo->dimensions[0];
	if (vtkInfo->operation == copyTo)
	{
		// Initialize the pointer to which it is to be copied to
		vtkInfo->dataPointer = (dataType **)malloc(sizeof(dataType *) * height);
		for (i = 0; i < height; i++)
		{
			vtkInfo->dataPointer[i] = (dataType *)malloc(sizeof(dataType) * (length * width));
		}
	}
	for (k = 0; k < height; k++)
	{
		for (i = 0; i < length; i++)
		{
			for (j = 0; j < width; j++)
			{
				xd = x_new(i, j, length);
				xyd = x_flat(i, j, k, length, width);
				if (vtkInfo->operation == copyFrom)
				{
					ptr[xyd] = static_cast<int>(vtkInfo->dataPointer[k][xd]);
				}
				else if (vtkInfo->operation == copyTo)
				{
					vtkInfo->dataPointer[k][xd] = static_cast<dataType>(ptr[xyd]);
				}
			}
		}
	}
}
//==============================================================================
void fillPtr(float * ptr, Vtk_File_Info * vtkInfo)
{
	size_t k, i, j, xd, xyd;
	int height = vtkInfo->dimensions[2], length = vtkInfo->dimensions[1], width = vtkInfo->dimensions[0];
	if (vtkInfo->operation == copyTo)
	{
		// Initialize the pointer to which it is to be copied to
		vtkInfo->dataPointer = (dataType **)malloc(sizeof(dataType *) * height);
		for (i = 0; i < height; i++)
		{
			vtkInfo->dataPointer[i] = (dataType *)malloc(sizeof(dataType) * (length * width));
		}
	}
	for (k = 0; k < height; k++)
	{
		for (i = 0; i < length; i++)
		{
			for (j = 0; j < width; j++)
			{
				xd = x_new(i, j, length);
				xyd = x_flat(i, j, k, length, width);
				if (vtkInfo->operation == copyFrom)
				{
					ptr[xyd] = static_cast<float>(vtkInfo->dataPointer[k][xd]);
				}
				else if (vtkInfo->operation == copyTo)
				{
					vtkInfo->dataPointer[k][xd] = static_cast<dataType>(ptr[xyd]);
				}
			}
		}
	}
}
//==============================================================================
void fillPtr(unsigned int * ptr, Vtk_File_Info * vtkInfo)
{
	size_t k, i, j, xd, xyd;
	int height = vtkInfo->dimensions[2], length = vtkInfo->dimensions[1], width = vtkInfo->dimensions[0];
	if (vtkInfo->operation == copyTo)
	{
		// Initialize the pointer to which it is to be copied to
		vtkInfo->dataPointer = (dataType **)malloc(sizeof(dataType *) * height);
		for (i = 0; i < height; i++)
		{
			vtkInfo->dataPointer[i] = (dataType *)malloc(sizeof(dataType) * (length * width));
		}
	}
	for (k = 0; k < height; k++)
	{
		for (i = 0; i < length; i++)
		{
			for (j = 0; j < width; j++)
			{
				xd = x_new(i, j, length);
				xyd = x_flat(i, j, k, length, width);
				if (vtkInfo->operation == copyFrom)
				{
					ptr[xyd] = static_cast<unsigned int>(vtkInfo->dataPointer[k][xd]);
				}
				else if (vtkInfo->operation == copyTo)
				{
					vtkInfo->dataPointer[k][xd] = static_cast<dataType>(ptr[xyd]);
				}
			}
		}
	}
}
//==============================================================================
void fillPtr(unsigned char * ptr, Vtk_File_Info * vtkInfo)
{
	size_t k, i, j, xd, xyd;
	int height = vtkInfo->dimensions[2], length = vtkInfo->dimensions[1], width = vtkInfo->dimensions[0];
	if (vtkInfo->operation == copyTo)
	{
		// Initialize the pointer to which it is to be copied to
		vtkInfo->dataPointer = (dataType **)malloc(sizeof(dataType *) * height);
		for (i = 0; i < height; i++)
		{
			vtkInfo->dataPointer[i] = (dataType *)malloc(sizeof(dataType) * (length * width));
		}
	}
	for (k = 0; k < height; k++)
	{
		for (i = 0; i < length; i++)
		{
			for (j = 0; j < width; j++)
			{
				xd = x_new(i, j, length);
				xyd = x_flat(i, j, k, length, width);
				if (vtkInfo->operation == copyFrom)
				{
					ptr[xyd] = static_cast<unsigned char>(vtkInfo->dataPointer[k][xd]);
				}
				else if (vtkInfo->operation == copyTo)
				{
					vtkInfo->dataPointer[k][xd] = static_cast<dataType>(ptr[xyd]);
				}
			}
		}
	}
}
//==============================================================================