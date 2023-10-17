#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#ifndef SEGMENTATIO3D_ATLAS
#define SEGMENTATION3D_ATLAS
	//==============================================================================
	/*
	* Header file to contain general functions, variables, structs
	* These can be used in different segmentation algorithms
	* Avoids redefinitions elsewhere
	*/
	//==============================================================================
	// Includes
#include "segmentation3D_common.h"
//==============================================================================
// Macro's
//==============================================================================
// STRUCTs
// Image Container and Properties
typedef struct {
	char name[50];
	// Data Dimensions
	size_t height, length, width; // Absolute Dimension
	unsigned char** dta_v;

	dataType** dta_u, ** dta_u0, ** dta_uq, ** dta_tmp, ** dta_tmp2, ** dta_tmp3, ** dta_edgePointer, ** priorShape, ** estShape; // Calculated Data Containers
	dataType** dta_h, ** dta_l, ** dta_w; // Calculated Central Data Containers - Vertical, Forward, Backward
	// Original Image Container
	dataType** dta;
	int obj_points; // Count contour points = 0.0
	bool est_start;
	dataType imgBackground, imgForeground;
} ABSContainer;
//==============================================================================
// Gradient Container and Properties
typedef struct {
	// Calculated data Gradients
	dataType** ae, ** aw, ** an, ** as, ** at, ** ab;
	// Original data gradients
	dataType** ge, ** gw, ** gn, ** gs, ** gt, ** gb, ** gd;
	// Holder of 1 + gradients - ap
	dataType** ap;
} GradData;
//==============================================================================
// Container Ap, Bp, pq for the Atlas segment and Properties
typedef struct {
	// Holder of bp
	dataType** bp_e, ** bp_w, ** bp_n, ** bp_s, ** bp_b, ** bp_t;
	// Holder of cp, dp
	dataType** cp, ** dp;
	// Original Gradient Holders g2
	dataType** gd;
} AtlasData;
//==============================================================================
// Segmentation parameters 
typedef struct {
	size_t maxNoGaussIterations, reflectionLength, storeFileFrequency, timeSteps;
	dataType episonRegularization, tau, edgeDetectorCoef, inflationCoeff, toleranceGausedel, toleranceSegmentation, omegaGaussedel;
	dataType lamda, eta, eta; // control parameters for expanding and edge term influences respectively
	Point3D hSpacing; // pixel spacing - x,y,z
	dataType mu1, mu2; // control parameters for the external force and the curvature terms
	dataType d1, d2; // distance used to determing what gamma value to set before using the atlas
} Segmentation_Paramereters;
//==============================================================================
// Struct for mean shape, eigenvectors, eigenvalues, princomp
typedef struct {
	// Mean Shape, eigenvector
	dataType** dtaMean, ** eigenvectors;
	// Distance maps
	dataType** dtaMeanDist;
	// Principal component
	dataType* eigenvalues;
	int princomp; enum EstimationMethod estMethod;
} PCAData;
//==============================================================================
// ENUMs
enum EstimationMethod { MIN_Energy = 1, MIN_Probability };
//==============================================================================
// Global Variables
//==============================================================================
// Function Prototypes
void segmentation3D_Ap_coef(ABSContainer* dta3D, GradData* grad3D, AtlasData* atls3D, size_t height, size_t length, size_t width, T** priorShape, T lambda, T eta, T zeta, T h, size_t p, T epsilon, T tau);
//==============================================================================
void gmcf3D_atlas(ABSContainer* dta3D, GradData* grad3D, AtlasData* atls3D, size_t height, size_t length, size_t width, size_t p);
//==============================================================================
// Atlas Segmentation Model interface function
void atlasSegmentationModel(Segmentation_Paramereters segParameters, size_t imageHeight, size_t imageLength, size_t imageWidth);
//==============================================================================
#endif // !SEGMENTATIO3D_COMMON

#ifdef __cplusplus
}
#endif