
#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#ifndef COMMON_FUNCTIONS
#define COMMON_FUNCTIONS
	//==============================================================================
	//debug constants
#ifndef MEASURE_TIME
#define MEASURE_TIME
#endif
// Show/Hide Output in Console
#ifndef CONSOLE_OUTPUT
#define CONSOLE_OUTPUT
#endif
#undef CONSOLE_OUTPUT
//==============================================================================
// Typedefs
// Data Type

#ifndef USE_DOUBLE_PRECISION
	typedef float dataType;
#else
	typedef double dataType;
#endif // USE_DOUBLE_PRECISION

	//==============================================================================
	// Includes
#include <stddef.h>
#include <omp.h>
#include <stdbool.h>
//==============================================================================
// MACROs
//==============================================================================
// STRUCTs

// Common 2D Points - {x,y}
	typedef struct ptstruct{
		dataType x;
		dataType y;
	} Point2D;

	// 2D point extended by flag indicating, if the point is end point (1st or last)
	typedef struct {
		union {
			struct ptstruct;
			Point2D pt;
		};
		bool isEndPoint;
		void * prev_point;
		void * next_point;
	} CurvePoint2D;

	// 2D Curve - list of CurvePoint2D
	typedef struct {
		CurvePoint2D* pPoints;
		size_t numPoints;
	} Curve2D;

// Common 3D Points - {x,y,z}
	typedef struct {
		dataType x, y, z;
	} Point3D;

	//Structure to handle image spacing
	typedef struct {
		dataType sx, sy, sz;
	} VoxelSpacing;

	//Matrix for rotation
	typedef struct {
		Point3D v1, v2, v3;
	}OrientationMatrix;

	typedef struct {
		dataType hx;
		dataType hy;
	} FiniteVolumeSize2D;

	// Image Container and Properties
	typedef struct {
		// Image Dimensions
		size_t height, length, width; // Absolute Dimension
		dataType** imageDataPtr; // Image Data Containers
		Point3D origin; // image origin
		VoxelSpacing spacing; // distance between pixels and distance between slice //--> voxel dimension
		OrientationMatrix orientation;
	} Image_Data;

	// Generate Random Points
	typedef struct {
		size_t k, xd, p;
	}Random3dPoints;

	typedef struct {
		size_t xd, p;
	}Random2dPoints;

	typedef struct {
		Point2D v1, v2;
	}OrientationMatrix2D;

	typedef struct {
		dataType sx, sy;
	} PixelSpacing;

	typedef struct {
		// Image Dimensions
		size_t height, width; // Absolute Dimension
		dataType* imageDataPtr; // Image Data Containers
		Point2D origin; // image origin
		PixelSpacing spacing; // pixel size
		OrientationMatrix2D orientation;
	} Image_Data2D;

	typedef struct {
		dataType min_data, max_data, mean_data, sd_data;
	} Statistics;

	//==============================================================================
	// Shapes Container
	typedef struct {
		// Shapes
		dataType** shp;
		size_t num; // shape number
	} Shapes;
	//==============================================================================
	// Enumeration
	enum FiniteDifference
	{
		FINITE_FORWARD = 1, FINITE_BACKWARD, FINITE_CENTRAL
	};
	//==============================================================================
	// FUNCTION PROTOTYPES
	/*
	* Calculate to 1D representation from 2D Array
	* x from column and row indices
	*/
	size_t x_new(const size_t rowIndex, const size_t columnIndex, const size_t rowLength);
	//==============================================================================
	/* Flattens 3D array to 1D array
	*/
	size_t x_flat(const size_t rowIndex, const size_t columnIndex, const size_t heightIndex, const size_t rowLength, const size_t columnLength);
	//==============================================================================
	/*
	* Reflection Function Zuska
	* toReflectImage is the data to perform reflection on
	* imageHeight is the actual Z dimension
	* imageLength is the actual X dimension
	* imageWidth is the actual Y dimension
	* P is the boundary constant/distance
	* Reflection does not include the current point but only its adjacent neighbor
	* Preserves Neumann Boundary condition - Target is for those methods
	*/
	void reflection3DB(dataType** toReflectImage, size_t imageHeight, size_t imageLength, size_t imageWidth, size_t p);
	//==============================================================================
	/*
	* Reflection Function Jozef
	* Reflection includes the current point together with its adjacent neighbor
	* Targets Distance Function types
	*/
	void reflection3D(dataType** toReflectImage, size_t imageHeight, size_t imageLength, size_t imageWidth);
	//==============================================================================
	/*
	* Edge Detector Calculation function
	*/
	dataType edgeDetector(dataType value, dataType coef);
	//==============================================================================
	/*
	void copyDataToExtendedArea(const dataType ** originalDataPtr, dataType ** extendedDataPtr, const size_t originalHeight, const size_t originalLength, const size_t originalWidth);

	The function copies data from source given by originalDataPtr
	representing data of dimensions given by originalHeight, originalLength and originalWidth
	to an array extendedDataPtr representing enlarged data (by 1 vx in each direction)
	of dimensions originalHeight + 2, originalLength + 2 and originalWidth + 2
	*/
	void copyDataToExtendedArea(dataType** originalDataPtr, dataType** extendedDataPtr, const size_t originalHeight, const size_t originalLength, const size_t originalWidth);
	//==============================================================================
	void copyDataToReducedArea(dataType** originalDataPtr, const dataType** extendedDataPtr, const size_t originalHeight, const size_t originalLength, const size_t originalWidth);
	//==============================================================================
	// Copy from one pointer to another pointer
	void copyDataToAnotherArray(dataType** source, dataType** destination, size_t height, size_t length, size_t width);
	//==============================================================================
	void rescaleNewRange(dataType** imageDataPtr, size_t imageLength, size_t imageWidth, size_t imageHeight, dataType minNew, dataType maxNew, dataType max_dta, dataType min_dta);
	//==============================================================================
	typedef struct {
		size_t k_min, i_min, j_min, k_max, i_max, j_max;
	} ClipBox;
	//==============================================================================
	// Calc. centroid of image data
	void centroidImage(dataType** imageDataPtr, dataType* centroid, size_t imageHeight, size_t imageLength, size_t imageWidth, dataType imageBackground);
	void centroidClipBox(dataType* centroid, ClipBox coord, dataType** imageDataPtr, size_t imageLength, dataType imageBackground);
	//==============================================================================

	//==============================================================================
	void copyDataToAnother2dArray(dataType* source, dataType* destination, size_t imageHeight, size_t imageWidth);
	//==============================================================================
	void copyDataTo2dExtendedArea(dataType* originalDataPtr, dataType* extendedDataPtr, const size_t originalHeight, const size_t originalWidth);
	//==============================================================================
	void copyDataTo2dReducedArea(dataType* originalDataPtr, const dataType* extendedDataPtr, const size_t originalHeight, const size_t originalWidth);
	//==============================================================================
	void reflection2D(dataType* toReflectImage, size_t imageHeight, size_t imageWidth);
	//==============================================================================
	double getPoint2DDistance(const Point2D a, const Point2D b);
	
	/// <summary>
	/// The points of pCurve are copied to pArray. The function expects same length of the particular objects (curve and array)
	/// </summary>
	/// <param name="pCurve">Pointer to the curve to be copied</param>
	/// <param name="array">Alocated array of expected size length (same as the curve) * 2</param>
	/// <returns></returns>
	bool copyCurve2DPointsToArray(const Curve2D * pCurve, dataType ** pArray);

	/// <summary>
	/// The function check, if the given curve is closed or not
	/// </summary>
	/// <param name="pcurve">Input curve to be checked</param>
	/// <returns>Returns true, if the curve is closed, otherwise false.</returns>
	bool isCurveClosed(const Curve2D * pcurve);

	/// <summary>
	/// Calculates gradient in 2D point given by central difference on input data
	/// </summary>
	/// <param name="image_data">Input image data</param>
	/// <param name="ind_x">x coordinate of the finite volume to calculate the gradient component</param>
	/// <param name="ind_y">y coordinate of the finite volume to calculate the gradient component</param>
	/// <param name="sz">size of finite volumes</param>
	/// <param name="grad">output - calculated gradient</param>
	/// <returns>True, if it was possible to estimate gradient</returns>
	bool getGradient2D(Image_Data2D image_data, const size_t ind_x, const size_t ind_y, const FiniteVolumeSize2D sz, Point2D * grad);

	/// <summary>
	/// The function returns the distance to given point from orgin (0,0) - in other words, calculated a norm of the given vector 
	/// </summary>
	/// <param name="pt">Given input point</param>
	/// <returns>Returns the result of (sqrt(pt.x * pt.x + pt.y * pt.y))</returns>
	dataType norm(const Point2D pt);

	//==============================================================================
	double getPoint3DDistance(const Point3D a, const Point3D b);
	//==============================================================================
	/*
	* Point3D getPointWithTheHighestValue(dataType** distanceMapPtr, const size_t length, const size_t width, const size_t height)
	* distanceMapPtr : pointer contaning the computed distance for each pixel
	* lenght, width, height : image dimension
	* The function return the coordinates of the voxel with the higest value
	*/
	Point3D getPointWithTheHighestValue(dataType** distanceMapPtr, const size_t length, const size_t width, const size_t height);
#endif // !COMMON_FUNCTIONS

#ifdef __cplusplus
}
#endif