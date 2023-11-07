#include "segmentation3D_atlas.h"
#include "filter_params.h"
#include <limits.h>
//==============================================================================
// Function to calclate the gamma parameter - eq 91
dataType GammaFun(dataType value, dataType dist1, dataType dist2);
//==============================================================================
// Intitialize the atlas data containers
void initializeAtlasData(AtlasData* atls3D, size_t height, size_t length, size_t width, size_t reflexLength);
// Initialize Gradient data containers
void initializeGradData(GradData* grad3D, size_t imageHeight, size_t imageLength, size_t imageWidth, size_t reflexLength);
// Initialize the ABSContainer
void initializeABSContainer(ABSContainer* dta3D, size_t imageHeight, size_t imageLength, size_t imageWidth, size_t reflexLength);
// Estimate shape function by minimization of the energy
void estimatedShapeByEnergy(ABSContainer* dta3D, dataType** dtaMeanShape, dataType** shape, dataType** currentSeg, dataType* eigenvalues, dataType** eigenvectors, int princomp, size_t height, size_t length, size_t width, Affine_Parameter finalResults, dataType bound, dataType eps, int est_iterations, bool parallelize, int NUMTHREADS, dataType h, dataType tolerance);
// Estimate shape function by probability function
void estimateShapeByProbability(ABSContainer* dta3D, dataType** dtaMeanShape, dataType** shape, dataType** currentSeg, dataType* eigenvalues, dataType** eigenvectors, int princomp, size_t height, size_t length, size_t width, Affine_Parameter finalResults, dataType bound, dataType eps, int est_iterations, bool parallelize, int NUMTHREADS, dataType h, dataType tolerance);
//==============================================================================
// Free Containers
void freeAtlasABSContainerGradData(ABSContainer* dta3D, AtlasData* atls3D, GradData* grad3D, size_t imageHeight, size_t reflexLength);
//==============================================================================
// Calculate the variance fun
void calcVariance(dataType** aPtr, dataType** bPtr, dataType distWithin, size_t height, size_t length, size_t width, size_t p, int step, dataType rho_zero_step);
//==============================================================================
// Calculate the rho value
dataType calcRho(dataType** aPtr, dataType** bPtr, dataType distWithin, size_t height, size_t length, size_t width, size_t p);
//==============================================================================
// 
dataType error_3D(dataType** aPtr, dataType** bPtr, size_t height, size_t length, size_t width, dataType h_val, size_t step);
//==============================================================================
// Curve functions for estimation
void extend_array(dataType** ext_array, dataType* in_array, int length, int w_size);
void moving_average(dataType** results, dataType* dtarray, int w_size, int length);
void est_curve(dataType** est_y, dataType* arr, dataType* d_arr, int count_b, int count_e, int* est_x, int length, int offset);
dataType curve_area(dataType* arr, int length, int interval);
int half_area(dataType* arr, int count_a, int length);
int find_min(dataType* arr, dataType* minimum, int begin, int length);
void diff_s(dataType** rlts, dataType* arr, int length);
//==============================================================================
/* Function implementations*/
//==============================================================================
dataType GammaFun(dataType value, dataType dist1, dataType dist2)
{
	if (fabs(value) > dist1)
	{
		return 0.;
	}
	else if (fabs(value) < dist2)
	{
		return 1.;
	}
	else
	{
		return 0.5;
	}
}
//==============================================================================
void initializeAtlasData(AtlasData* atls3D, size_t imageHeight, size_t imageLength, size_t imageWidth, size_t reflexLength)
{
	size_t lengthExt = imageLength + 2 * reflexLength + 1, widthExt = imageWidth + 2 * reflexLength + 1, heightExt = imageHeight + 2 * reflexLength + 1, dim2DExt = lengthExt * widthExt, i, j;
	//==============================================================================
	// Initial value
	dataType inititialValue = 0.0f;
	//==============================================================================
	atls3D = (AtlasData*)malloc(sizeof(AtlasData));
	//==============================================================================
	// Allocate memory
	// Bp
	(*atls3D).bp_e = (dataType**)malloc(heightExt * sizeof(dataType*));
	(*atls3D).bp_w = (dataType**)malloc(heightExt * sizeof(dataType*));
	(*atls3D).bp_s = (dataType**)malloc(heightExt * sizeof(dataType*));
	(*atls3D).bp_n = (dataType**)malloc(heightExt * sizeof(dataType*));
	(*atls3D).bp_b = (dataType**)malloc(heightExt * sizeof(dataType*));
	(*atls3D).bp_t = (dataType**)malloc(heightExt * sizeof(dataType*));
	// Cp. dp
	(*atls3D).cp = (dataType**)malloc(heightExt * sizeof(dataType*));
	(*atls3D).dp = (dataType**)malloc(heightExt * sizeof(dataType*));
	// G2
	(*atls3D).gd = (dataType**)malloc(heightExt * sizeof(dataType*));
	//==============================================================================
	// Column memory allocation
	for (i = 0; i < heightExt; i++)
	{
		// Bp
		(*atls3D).bp_e[i] = (dataType*)malloc(dim2DExt * sizeof(dataType));
		(*atls3D).bp_w[i] = (dataType*)malloc(dim2DExt * sizeof(dataType));
		(*atls3D).bp_s[i] = (dataType*)malloc(dim2DExt * sizeof(dataType));
		(*atls3D).bp_n[i] = (dataType*)malloc(dim2DExt * sizeof(dataType));
		(*atls3D).bp_b[i] = (dataType*)malloc(dim2DExt * sizeof(dataType));
		(*atls3D).bp_t[i] = (dataType*)malloc(dim2DExt * sizeof(dataType));
		// Cp, dp
		(*atls3D).cp[i] = (dataType*)malloc(dim2DExt * sizeof(dataType));
		(*atls3D).dp[i] = (dataType*)malloc(dim2DExt * sizeof(dataType));
		// G2
		(*atls3D).gd[i] = (dataType*)malloc(dim2DExt * sizeof(dataType));

		for (j = 0; j < dim2DExt; j++)
		{
			// Bp
			(*atls3D).bp_e[i][j] = inititialValue;
			(*atls3D).bp_w[i][j] = inititialValue;
			(*atls3D).bp_s[i][j] = inititialValue;
			(*atls3D).bp_n[i][j] = inititialValue;
			(*atls3D).bp_b[i][j] = inititialValue;
			(*atls3D).bp_t[i][j] = inititialValue;
			// Cp
			(*atls3D).cp[i][j] = inititialValue;
			(*atls3D).dp[i][j] = inititialValue;
			// G2
			(*atls3D).gd[i][j] = inititialValue;
		}
	}
	//==============================================================================

}
//==============================================================================
void initializeGradData(GradData* grad3D, size_t imageHeight, size_t imageLength, size_t imageWidth, size_t reflexLength)
{
	size_t lengthExt = imageLength + 2 * reflexLength + 1, widthExt = imageWidth + 2 * reflexLength + 1, heightExt = imageHeight + 2 * reflexLength + 1, dim2DExt = lengthExt * widthExt, i, j;
	//==============================================================================
	// Initial value
	dataType inititialValue = 0.0f;
	//==============================================================================
	grad3D = (GradData*)malloc(sizeof(GradData));
	//==============================================================================
	// Allocate memory
	// Containers for tau/m(p) * g1|\nable U_p | * (1 / {|\nabla U_pq|)
	(*grad3D).ae = (dataType**)malloc(heightExt * sizeof(dataType*));
	(*grad3D).aw = (dataType**)malloc(heightExt * sizeof(dataType*));
	(*grad3D).an = (dataType**)malloc(heightExt * sizeof(dataType*));
	(*grad3D).as = (dataType**)malloc(heightExt * sizeof(dataType*));
	(*grad3D).at = (dataType**)malloc(heightExt * sizeof(dataType*));
	(*grad3D).ab = (dataType**)malloc(heightExt * sizeof(dataType*));
	// Containers for g1|\nable U_p |
	(*grad3D).ge = (dataType**)malloc(heightExt * sizeof(dataType*));
	(*grad3D).gw = (dataType**)malloc(heightExt * sizeof(dataType*));
	(*grad3D).gn = (dataType**)malloc(heightExt * sizeof(dataType*));
	(*grad3D).gs = (dataType**)malloc(heightExt * sizeof(dataType*));
	(*grad3D).gt = (dataType**)malloc(heightExt * sizeof(dataType*));
	(*grad3D).gb = (dataType**)malloc(heightExt * sizeof(dataType*));

	(*grad3D).gd = (dataType**)malloc(heightExt * sizeof(dataType*));
	//==============================================================================
	// Column memory allocation
	for (i = 0; i < heightExt; i++)
	{
		(*grad3D).ae[i] = (dataType*)malloc(dim2DExt * sizeof(dataType));
		(*grad3D).aw[i] = (dataType*)malloc(dim2DExt * sizeof(dataType));
		(*grad3D).an[i] = (dataType*)malloc(dim2DExt * sizeof(dataType));
		(*grad3D).as[i] = (dataType*)malloc(dim2DExt * sizeof(dataType));
		(*grad3D).at[i] = (dataType*)malloc(dim2DExt * sizeof(dataType));
		(*grad3D).ab[i] = (dataType*)malloc(dim2DExt * sizeof(dataType));

		(*grad3D).ge[i] = (dataType*)malloc(dim2DExt * sizeof(dataType));
		(*grad3D).gw[i] = (dataType*)malloc(dim2DExt * sizeof(dataType));
		(*grad3D).gn[i] = (dataType*)malloc(dim2DExt * sizeof(dataType));
		(*grad3D).gs[i] = (dataType*)malloc(dim2DExt * sizeof(dataType));
		(*grad3D).gt[i] = (dataType*)malloc(dim2DExt * sizeof(dataType));
		(*grad3D).gb[i] = (dataType*)malloc(dim2DExt * sizeof(dataType));

		(*grad3D).gd[i] = (dataType*)malloc(dim2DExt * sizeof(dataType));

		for (j = 0; j < dim2DExt; j++)
		{
			(*grad3D).ae[i][j] = inititialValue;
			(*grad3D).aw[i][j] = inititialValue;
			(*grad3D).an[i][j] = inititialValue;
			(*grad3D).as[i][j] = inititialValue;
			(*grad3D).at[i][j] = inititialValue;
			(*grad3D).ab[i][j] = inititialValue;

			(*grad3D).ge[i][j] = inititialValue;
			(*grad3D).gw[i][j] = inititialValue;
			(*grad3D).gn[i][j] = inititialValue;
			(*grad3D).gs[i][j] = inititialValue;
			(*grad3D).gt[i][j] = inititialValue;
			(*grad3D).gb[i][j] = inititialValue;

			(*grad3D).gd[i][j] = inititialValue;
		}
	}
	//==============================================================================

}
//==============================================================================
// Initialize the ABSContainer
void initializeABSContainer(ABSContainer* dta3D, size_t imageHeight, size_t imageLength, size_t imageWidth, size_t reflexLength)
{
	size_t lengthExt = imageLength + 2 * reflexLength + 1, widthExt = imageWidth + 2 * reflexLength + 1, heightExt = imageHeight + 2 * reflexLength + 1, dim2DExt = lengthExt * widthExt, i, j;
	//==============================================================================
	// Initial value
	dataType inititialValue = 0.0f;
	//==============================================================================
	dta3D = (ABSContainer*)malloc(sizeof(ABSContainer));
	//==============================================================================
	// Assisgn the dimensions
	(*dta3D).height = imageHeight;
	(*dta3D).length = imageLength;
	(*dta3D).width = imageWidth;
	//==============================================================================
	// Center Containers
	(*dta3D).dta_h = (dataType**)malloc(heightExt * sizeof(dataType*));
	(*dta3D).dta_l = (dataType**)malloc(heightExt * sizeof(dataType*));
	(*dta3D).dta_w = (dataType**)malloc(heightExt * sizeof(dataType*));
	//==============================================================================
	// Data containers
	(*dta3D).dta_v = (unsigned char**)malloc(heightExt * sizeof(unsigned char*));
	(*dta3D).dta_u = (dataType**)malloc(heightExt * sizeof(dataType*));
	(*dta3D).dta = (dataType**)malloc(heightExt * sizeof(dataType*));
	(*dta3D).dta_u0 = (dataType**)malloc(heightExt * sizeof(dataType*));
	(*dta3D).dta_uq = (dataType**)malloc(heightExt * sizeof(dataType*));
	(*dta3D).dta_tmp = (dataType**)malloc(heightExt * sizeof(dataType*));
	(*dta3D).dta_tmp2 = (dataType**)malloc(heightExt * sizeof(dataType*));
	(*dta3D).dta_tmp3 = (dataType**)malloc(heightExt * sizeof(dataType*));
	(*dta3D).dta_edgePointer = (dataType**)malloc(heightExt * sizeof(dataType*));
	//==============================================================================
	(*dta3D).imgBackground = 1.0f, (*dta3D).imgForeground = 0.0f;
	//==============================================================================
	// Column memory allocations
	for (i = 0; i < heightExt; i++)
	{
		// Center Containers
		(*dta3D).dta_h[i] = (dataType*)malloc(dim2DExt * sizeof(dataType));
		(*dta3D).dta_l[i] = (dataType*)malloc(dim2DExt * sizeof(dataType));
		(*dta3D).dta_w[i] = (dataType*)malloc(dim2DExt * sizeof(dataType));
		// Data containers
		(*dta3D).dta_v[i] = (unsigned char*)malloc(heightExt * sizeof(unsigned char));
		(*dta3D).dta_u[i] = (dataType*)malloc(dim2DExt * sizeof(dataType));
		(*dta3D).dta[i] = (dataType*)malloc(dim2DExt * sizeof(dataType));
		(*dta3D).dta_u0[i] = (dataType*)malloc(dim2DExt * sizeof(dataType));
		(*dta3D).dta_uq[i] = (dataType*)malloc(dim2DExt * sizeof(dataType));
		(*dta3D).dta_tmp[i] = (dataType*)malloc(dim2DExt * sizeof(dataType));
		(*dta3D).dta_tmp2[i] = (dataType*)malloc(dim2DExt * sizeof(dataType));
		(*dta3D).dta_tmp3[i] = (dataType*)malloc(dim2DExt * sizeof(dataType));
		(*dta3D).dta_edgePointer[i] = (dataType*)malloc(dim2DExt * sizeof(dataType));

		for (j = 0; j < dim2DExt; j++)
		{
			// Center Containers
			(*dta3D).dta_h[i][j] = inititialValue;
			(*dta3D).dta_l[i][j] = inititialValue;
			(*dta3D).dta_w[i][j] = inititialValue;
			// Data containers
			(*dta3D).dta_v[i][j] = 0;
			(*dta3D).dta_u[i][j] = inititialValue;
			(*dta3D).dta[i][j] = inititialValue;
			(*dta3D).dta_u0[i][j] = inititialValue;
			(*dta3D).dta_uq[i][j] = inititialValue;
			(*dta3D).dta_tmp[i][j] = inititialValue;
			(*dta3D).dta_tmp2[i][j] = inititialValue;
			(*dta3D).dta_tmp3[i][j] = inititialValue;
			(*dta3D).dta_edgePointer[i][j] = inititialValue;
		}
	}
	//==============================================================================
}
//==============================================================================
// Estimate shape functions
//==============================================================================
void estimatedShapeByEnergy(ABSContainer* dta3D,dataType** dtaMeanShape, dataType** shape, dataType** currentSeg, dataType* eigenvalues, dataType** eigenvectors, int princomp, size_t height, size_t length, size_t width, Affine_Parameter finalResults, dataType bound, dataType eps, int est_iterations, bool parallelize, int NUMTHREADS, dataType h, dataType tolerance)
{
	int k, i, j, xd, xyd, l;
	//==============================================================================
	size_t dim3D = height * length * width, dim2D = length * width;
	//==============================================================================
	// Create k by 1 bounded eigen vectors
	dataType* bound_eigvalues = (dataType*)malloc(sizeof(dataType) * princomp); // K by 1
	for (i = 0; i < princomp; i++)
	{
		bound_eigvalues[i] = bound * sqrtf(fabsf(eigenvalues[i]));
	}
	//==============================================================================
	// Transpose the eigenvectors
	dataType** trans_eigvectors = (dataType**)malloc(sizeof(dataType*) * princomp); // K by D
	for (i = 0; i < princomp; i++)
	{
		trans_eigvectors[i] = (dataType*)malloc(sizeof(dataType) * dim3D);
	}
	transpose(trans_eigvectors, eigenvectors, dim3D, princomp);
	//==============================================================================
	// Difference btn shape and mean shape - D by 1 ---> Center the shape
	dataType* shp_diff = (dataType*)malloc(sizeof(dataType) * dim3D); // D by 1
	for (k = 0; k < height; k++)
	{
		for (i = 0; i < length; i++)
		{
			for (j = 0; j < width; j++)
			{
				// 2D to 1D representation for i, j
				xd = x_new(i, j, length);
				// 3D to 1D flattening
				xyd = x_flat(i, j, k, length, width);
				// Difference btweeen shape and mean shape
				shp_diff[xyd] = shape[k][xd] - dtaMeanShape[k][xd];
			}
		}
	}
	//==============================================================================
	// Projection of the shape - k by 1
	dataType* shp_projection = (dataType*)malloc(sizeof(dataType) * princomp); // k by 1
	// Multiplication
	for (i = 0; i < princomp; i++)
	{
		shp_projection[i] = 0;
	}
	//==============================================================================
	// Estimated shape from atlas
	for (l = 0; l < princomp; l++)
	{
		for (k = 0; k < height; k++)
		{
			for (i = 0; i < length; i++)
			{
				for (j = 0; j < width; j++)
				{
					// 2D to 1D representation for i, j
					xd = x_new(i, j, length);
					// 3D to 1D flattening
					xyd = x_flat(i, j, k, length, width);
					(*dta3D).estShape[k][xd] = (eigenvectors[xyd][l] * (shp_projection[l])) + dtaMeanShape[k][xd];
				}
			}
		}
	}
	//==============================================================================
	// Gradient of energy
	dataType* grad_energy = (dataType*)malloc(sizeof(dataType) * princomp); // k by 1
	for (l = 0; l < princomp; l++)
	{
		grad_energy[l] = 0;
		for (k = 0; k < height; k++)
		{
			for (i = 0; i < length; i++)
			{
				for (j = 0; j < width; j++)
				{
					// 2D to 1D representation for i, j
					xd = x_new(i, j, length);
					// 3D to 1D flattening
					xyd = x_flat(i, j, k, length, width);
					//grad_energy[l] += (2 * (trans_eigvectors[l][xyd] * ((*dta3D).priorShape[k][xd] - currentSeg[k][xd]))) / dim3D;
					// Difference btn est. and reg. shapes
					grad_energy[l] += (2 * (trans_eigvectors[l][xyd] * ((*dta3D).estShape[k][xd] - shape[k][xd]))) / dim3D;
				}
			}
		}
	}
	//==============================================================================
	// Minimization of Energy
	int iteration = 0, steps = est_iterations;
	dataType energy_param = errorCalc((*dta3D).estShape, shape, height, length, width, h);
	//==============================================================================
	printf("\nMEnergy minimized at begin is %.8f\n", energy_param);
	//==============================================================================
	// Obtain the minimum shape projection
	bool stp_condition = false;
	while (!stp_condition)
	{
		if (energy_param < tolerance || iteration == steps)
		{
			stp_condition = true;
		}
		else
		{
			// Estimated shape from atlas
			if (parallelize)
			{
				omp_set_dynamic(0);
				omp_set_num_threads(NUMTHREADS);
#pragma omp parallel
				{
#pragma omp for private(l) schedule(static) nowait
					for (l = 0; l < princomp; l++)
					{
						shp_projection[l] = shp_projection[l] - eps * grad_energy[l];
						// Boundary check
						if (shp_projection[l] > bound_eigvalues[l])
						{
							shp_projection[l] = bound_eigvalues[l];
						}
						else if (shp_projection[l] < (-1 * bound_eigvalues[l]))
						{
							shp_projection[l] = -1 * bound_eigvalues[l];
						}
					}
				}
			}
			else
			{
				for (l = 0; l < princomp; l++)
				{
					shp_projection[l] = shp_projection[l] - eps * grad_energy[l];
					// Boundary check
					if (shp_projection[l] > bound_eigvalues[l])
					{
						shp_projection[l] = bound_eigvalues[l];
					}
					else if (shp_projection[l] < (-1 * bound_eigvalues[l]))
					{
						shp_projection[l] = -1 * bound_eigvalues[l];
					}
				}
			}
			// Estimated shape from atlas
			for (l = 0; l < princomp; l++)
			{
				for (k = 0; k < height; k++)
				{
					for (i = 0; i < length; i++)
					{
						for (j = 0; j < width; j++)
						{
							// 2D to 1D representation for i, j
							xd = x_new(i, j, length);
							// 3D to 1D flattening
							xyd = x_flat(i, j, k, length, width);
							(*dta3D).estShape[k][xd] = (eigenvectors[xyd][l] * (shp_projection[l])) + dtaMeanShape[k][xd];
						}
					}
				}
			}
			//==============================================================================
			// Grad_energy
			// Estimated shape from atlas
			if (parallelize)
			{
				omp_set_dynamic(0);
				omp_set_num_threads(NUMTHREADS);
#pragma omp parallel
				{
#pragma omp for private(l,k,i,j,xd,xyd) schedule(static) nowait
					for (l = 0; l < princomp; l++)
					{
						grad_energy[l] = 0;
						for (k = 0; k < height; k++)
						{
							for (i = 0; i < length; i++)
							{
								for (j = 0; j < width; j++)
								{
									// 2D to 1D representation for i, j
									xd = x_new(i, j, length);
									// 3D to 1D flattening
									xyd = x_flat(i, j, k, length, width);
									//grad_energy[l] += (2 * (trans_eigvectors[l][xyd] * ((*dta3D).priorShape[k][xd] - currentSeg[k][xd]))) / dim3D;
									// Difference btn est. and reg. shapes
									grad_energy[l] += (2 * (trans_eigvectors[l][xyd] * ((*dta3D).estShape[k][xd] - shape[k][xd]))) / dim3D;
								}
							}
						}
					}
				}
			}
			else
			{
				for (l = 0; l < princomp; l++)
				{
					grad_energy[l] = 0;
					for (k = 0; k < height; k++)
					{
						for (i = 0; i < length; i++)
						{
							for (j = 0; j < width; j++)
							{
								// 2D to 1D representation for i, j
								xd = x_new(i, j, length);
								// 3D to 1D flattening
								xyd = x_flat(i, j, k, length, width);
								//grad_energy[l] += (2 * (trans_eigvectors[l][xyd] * ((*dta3D).priorShape[k][xd] - currentSeg[k][xd]))) / dim3D;
								// Difference btn est. and reg. shapes
								grad_energy[l] += (2 * (trans_eigvectors[l][xyd] * ((*dta3D).estShape[k][xd] - shape[k][xd]))) / dim3D;
							}
						}
					}
				}
			}
			//energy_param = errorCalc((*dta3D).priorShape, currentSeg, height, length, width, h);
			// Difference btn est. and reg. shapes
			energy_param = errorCalc((*dta3D).estShape, shape, height, length, width, h);
			//==============================================================================
			iteration++;
			//==============================================================================
		}
	}
	//==============================================================================
	printf("\nThe number of iterations to calc. prior shape is %d\n", iteration);
	//==============================================================================
	energy_param = errorCalc((*dta3D).estShape, shape, height, length, width, h);
	//==============================================================================
	printf("\nMEnergy minimized at end is %.8f\n", energy_param);
	// Copy Estimate shape to PriorShape (transformed shape)
	copyDataPointer((*dta3D).estShape, (*dta3D).priorShape, height, length, width);
	//==============================================================================
	// Transform the estimated shape
	scaleDown((*dta3D).priorShape, height, length, width, finalResults);
	//==============================================================================
	printf("\nShape Projection vector:\n");
	for (j = 0; j < princomp; j++) {
		printf("%.4f\n", shp_projection[j]);
	}
	//==============================================================================
	// Free Memory
	for (i = 0; i < princomp; i++)
	{
		free(trans_eigvectors[i]);
	}
	//free(trans_eigvectors);
	free(bound_eigvalues);
	free(shp_projection);
	free(grad_energy);
	free(shp_diff);
	shp_projection = NULL;
	trans_eigvectors = NULL;
	grad_energy = NULL;
	shp_diff = NULL;
	bound_eigvalues = NULL;
}
//==============================================================================
void estimateShapeByProbability(ABSContainer* dta3D, dataType** dtaMeanShape, dataType** shape, dataType** currentSeg, dataType* eigenvalues, dataType** eigenvectors, int princomp, size_t height, size_t length, size_t width, Affine_Parameter finalResults, dataType bound, dataType eps, int est_iterations, bool parallelize, int NUMTHREADS, dataType h, dataType tolerance)
{
	int k, i, j, xd, xyd;
	//==============================================================================
	size_t dim3D = height * length * width;
	//==============================================================================
	// Create a diagonal eigenvalues matrix inverse
	dataType** diag_eigvalues_inv = (dataType**)malloc(sizeof(dataType*) * princomp); // K by k
	for (i = 0; i < princomp; i++)
	{
		diag_eigvalues_inv[i] = (dataType*)malloc(sizeof(dataType) * princomp);
		for (j = 0; j < princomp; j++)
		{
			if (i == j)
			{
				diag_eigvalues_inv[i][j] = 1 / eigenvalues[j];
			}
			else
			{
				diag_eigvalues_inv[i][j] = 0;
			}
		}
	}
	//==============================================================================
	// Transpose the eigenvectors
	dataType** trans_eigvectors = (dataType**)malloc(sizeof(dataType*) * princomp); // K by D
	for (i = 0; i < princomp; i++)
	{
		trans_eigvectors[i] = (dataType*)malloc(sizeof(dataType) * dim3D);
	}
	transpose(trans_eigvectors, eigenvectors, dim3D, princomp);
	//==============================================================================
	// Difference btn shape and mean shape - D by 1
	dataType* shp_diff = (dataType*)malloc(sizeof(dataType) * dim3D); // D by 1
	for (k = 0; k < height; k++)
	{
		for (i = 0; i < length; i++)
		{
			for (j = 0; j < width; j++)
			{
				// 2D to 1D representation for i, j
				xd = x_new(i, j, length);
				// 3D to 1D flattening
				xyd = x_flat(i, j, k, length, width);
				// Difference btweeen shape and mean shape
				shp_diff[xyd] = shape[k][xd] - dtaMeanShape[k][xd];
			}
		}
	}
	//==============================================================================
	// Projection of the shape - k by 1
	dataType* shp_projection = (dataType*)malloc(sizeof(dataType) * princomp); // k by 1
	dataType* tmp_results = (dataType*)malloc(sizeof(dataType) * princomp); // k by 1
	//==============================================================================
	// Multiplication
	for (i = 0; i < princomp; i++)
	{
		shp_projection[i] = 0;
		for (j = 0; j < dim3D; j++)
		{
			shp_projection[i] += trans_eigvectors[i][j] * shp_diff[j];
		}
		if (shp_projection[i] > bound * sqrtf(fabsf(eigenvalues[i])))
		{
			shp_projection[i] = bound * sqrtf(fabsf(eigenvalues[i]));
		}
		else if (shp_projection[i] < -1 * bound * sqrtf(fabsf(eigenvalues[i])))
		{
			shp_projection[i] = -1 * bound * sqrtf(fabsf(eigenvalues[i]));
		}
	}
	//==============================================================================
	// Estimated shape parameter - 1 by 1
	dataType energy_param = 0;
	// Shape projection transpose
	dataType** shp_projection_trans = (dataType*)malloc(sizeof(dataType) * 1); // 1 by k
	for (i = 0; i < 1; i++)
	{
		shp_projection_trans[i] = (dataType*)malloc(sizeof(dataType) * princomp);
	}
	for (i = 0; i < princomp; i++)
	{
		shp_projection_trans[0][i] = shp_projection[i];
	}
	//==============================================================================
	for (i = 0; i < princomp; i++)
	{
		tmp_results[i] = 0.;
		for (j = 0; j < princomp; j++)
		{
			tmp_results[i] += diag_eigvalues_inv[i][j] * shp_projection[j];
		}
	}
	//==============================================================================
	for (i = 0; i < princomp; i++)
	{
		energy_param += shp_projection_trans[0][i] * tmp_results[i];
	}
	//==============================================================================
	printf("\nMinimized EShape at begin: %.8lf\n", energy_param);
	//==============================================================================
	bool stpCond = false;
	//==============================================================================
	// Minimization of Energy
	int iteration = 0, steps = est_iterations;
	//==============================================================================
	// Initialize pointers#
	dataType* est_shape = (dataType*)malloc(sizeof(dataType) * dim3D); // S - D by 1 components
	//==============================================================================
	dataType* save_results = (dataType*)malloc(sizeof(dataType) * (1 + princomp));
	//==============================================================================
	// Calc. Estimated U_k*a_k
	for (i = 0; i < dim3D; i++)
	{
		est_shape[i] = 0.;
		for (j = 0; j < princomp; j++)
		{
			// Estimated shape
			est_shape[i] += shp_projection[j] * eigenvectors[i][j];
		}
	}
	//==============================================================================
	// Construct shape prior model
	for (k = 0; k < height; k++)
	{
		for (i = 0; i < length; i++)
		{
			for (j = 0; j < width; j++)
			{
				// 2D to 1D representation for i, j
				xd = x_new(i, j, length);
				// 3D to 1D flattening
				xyd = x_flat(i, j, k, length, width);
				// Construct prior shape
				(*dta3D).estShape[k][xd] = dtaMeanShape[k][xd] + est_shape[xyd];
			}
		}
	}
	//==============================================================================
	while (!stpCond)
	{
		if (energy_param < tolerance || iteration == steps)
		{
			stpCond = true;
		}
		else
		{
			//==============================================================================
			if (parallelize)
			{
				omp_set_dynamic(0);
				omp_set_num_threads(NUMTHREADS);
#pragma omp parallel
				{
#pragma omp for private(i) schedule(static) nowait
					for (i = 0; i < princomp; i++)
					{
						shp_projection[i] = shp_projection[i] - eps * tmp_results[i];
						if (shp_projection[i] > bound * sqrtf(fabsf(eigenvalues[i])))
						{
							shp_projection[i] = bound * sqrtf(fabsf(eigenvalues[i]));
						}
						else if (shp_projection[i] < -1 * bound * sqrtf(fabsf(eigenvalues[i])))
						{
							shp_projection[i] = -1 * bound * sqrtf(fabsf(eigenvalues[i]));
						}
					}
				}
			}
			else // Sequential
			{
				for (i = 0; i < princomp; i++)
				{
					shp_projection[i] = shp_projection[i] - eps * tmp_results[i];
					if (shp_projection[i] > bound * sqrtf(fabsf(eigenvalues[i])))
					{
						shp_projection[i] = bound * sqrtf(fabsf(eigenvalues[i]));
					}
					else if (shp_projection[i] < -1 * bound * sqrtf(fabsf(eigenvalues[i])))
					{
						shp_projection[i] = -1 * bound * sqrtf(fabsf(eigenvalues[i]));
					}
				}
			}
			//==============================================================================
			// Update Shape Projection T
			for (i = 0; i < princomp; i++)
			{
				shp_projection_trans[0][i] = shp_projection[i];
			}
			// Update the tmp results
			for (i = 0; i < princomp; i++)
			{
				tmp_results[i] = 0.;
				for (j = 0; j < princomp; j++)
				{
					tmp_results[i] += diag_eigvalues_inv[i][j] * shp_projection[j];
				}
			}
			//==============================================================================
			// Error
			energy_param = 0;
			for (i = 0; i < princomp; i++)
			{
				energy_param += shp_projection_trans[0][i] * tmp_results[i];
			}
			//==============================================================================
			save_results[0] = energy_param;
			for (i = 0; i < princomp; i++)
			{
				save_results[i + 1] = shp_projection[i];
			}
			//==============================================================================
			// Calc. Estimated U_k*a_k
			for (i = 0; i < dim3D; i++)
			{
				est_shape[i] = 0.;
				for (j = 0; j < princomp; j++)
				{
					// Estimated shape
					est_shape[i] += shp_projection[j] * eigenvectors[i][j];
				}
			}
			//==============================================================================
			// Construct shape prior model
			if (parallelize)
			{
				omp_set_dynamic(0);
				omp_set_num_threads(NUMTHREADS);
#pragma omp parallel
				{
#pragma omp for private(k,i,j,xd,xyd) schedule(static) nowait
					for (k = 0; k < height; k++)
					{
						for (i = 0; i < length; i++)
						{
							for (j = 0; j < width; j++)
							{
								// 2D to 1D representation for i, j
								xd = x_new(i, j, length);
								// 3D to 1D flattening
								xyd = x_flat(i, j, k, length, width);
								// Construct prior shape
								(*dta3D).estShape[k][xd] = dtaMeanShape[k][xd] + est_shape[xyd];
							}
						}
					}
				}
			}
			else // Sequential
			{
				for (k = 0; k < height; k++)
				{
					for (i = 0; i < length; i++)
					{
						for (j = 0; j < width; j++)
						{
							// 2D to 1D representation for i, j
							xd = x_new(i, j, length);
							// 3D to 1D flattening
							xyd = x_flat(i, j, k, length, width);
							// Construct prior shape
							(*dta3D).estShape[k][xd] = dtaMeanShape[k][xd] + est_shape[xyd];
						}
					}
				}
			}
			//==============================================================================
			iteration++;
			//==============================================================================
		}
	}
	//==============================================================================
	printf("\nThe number of iterations to calc. prior shape is %d\n", iteration);
	printf("\nMinimized EShape at end is: %.8lf\n", energy_param);
	//==============================================================================
	// Copy Estimate shape to PriorShape (transformed shape)
	copyDataPointer((*dta3D).estShape, (*dta3D).priorShape, height, length, width);
	//==============================================================================
	// Transform the estimated shape
	scaleDown((*dta3D).priorShape, height, length, width, finalResults);
	//==============================================================================
	printf("Shape Projection vector:\n");
	for (j = 0; j < princomp; j++) {
		printf("%.6f\n", shp_projection[j]);
	}
	//==============================================================================
	free(est_shape);
	est_shape = NULL;
	//==============================================================================
	free(save_results);
	save_results = NULL;
	//==============================================================================
	// Free Memory
	for (i = 0; i < princomp; i++)
	{
		free(diag_eigvalues_inv[i]);
		free(trans_eigvectors[i]);
	}
	for (i = 0; i < 1; i++)
	{
		free(shp_projection_trans[i]);
	}
	free(diag_eigvalues_inv);
	free(shp_projection_trans);
	free(trans_eigvectors);
	free(tmp_results);
	free(shp_projection);
	free(shp_diff);

	tmp_results = NULL;
	shp_projection = NULL;
	diag_eigvalues_inv = NULL;
	trans_eigvectors = NULL;
	shp_diff = NULL;
	shp_projection_trans = NULL;
	//==============================================================================
}
//==============================================================================
// Estimate the current segmentation from the atlas function
dataType estmateSegmentation(ABSContainer* dta3D, PCAData* pcaParam, dataType** currentSegmentation, size_t height, size_t length, size_t width, size_t reflexLength, dataType foreground, dataType background, dataType tol_s, dataType h, Registration_Params regParams, Optimization_Method optMethod)
{
	//==============================================================================
	size_t dim2D = length + 2 * reflexLength + 1 * width + 2 * reflexLength + 1;
	//==============================================================================
	// Params
	size_t k, smooth_level = 300;
	dataType error = INT_MAX, tol_est = 0.05 /*0.12, 0.15*/, step_size = 0.25;
	char st[160];
	for (k = 0; k < 160; k++)
	{
		st[k] = '=';
	}
	st[159] = '\0';
	//==============================================================================
	printf("%s\n", st);
	printf("Begin estimation of the current segmentation shape\n");
	printf("%s\n\n", st);
	// Check for Null struct before - only begin estimation when segmentation has slowed down
	if (!(*dta3D).est_start)
	{
		printf("Exiting estimation of current segmentation, NULL pca parameters\n");
		goto ExitCond;
	}
	else
	{
		// Register the current segmentation
		dataType** regCurrentSegmentation = (dataType**)malloc((height + 2 * reflexLength + 1) * sizeof(dataType*)); // Stores the registered segmentation
		for (k = 0; k < height + 2 * reflexLength + 1; k++)
		{
			regCurrentSegmentation[k] = (dataType*)malloc(dim2D * sizeof(dataType));
		}
		// Initialize to zeros
		initialize3dArrayD(regCurrentSegmentation, length + 2 * reflexLength + 1, width + 2 * reflexLength + 1, height + 2 * reflexLength + 1, foreground);
		//==============================================================================
		// Copy current segmentation to tmp pointer
		copyDataPointer(currentSegmentation, (*dta3D).dta_tmp3, height, length, width);
		//==============================================================================
		// Create a bin shape from current segmentation dist. map - used in registration
		double threshold = 0.7;
		BinaryMap((*dta3D).dta_tmp3, height, length, width, threshold, foreground, background);
		//==============================================================================
		// Apply smoothing to the bin map of current segmentation - level same as applied to the Mean shape when it was created
		dataType stepsize = 2.0, time_steps = (size_t)smooth_level * 1.5, tolerance = tol_s;
		smoothData((*dta3D).dta_tmp3, height, length, width, reflexLength, stepsize, time_steps, tolerance, 1);
		// Rescale 0 - 1
		rescaleDta((*dta3D).dta_tmp3, height, length, width, 0, 0);
		// Registration - Inputs are distance Maps for current segmentation and mean shape
		run_registration((*pcaParam).dtaMean, (*dta3D).dta_tmp3, regCurrentSegmentation, height, length, width, regParams, optMethod);
		//==============================================================================
		reflection3D(regCurrentSegmentation, height - reflexLength, length - reflexLength, width - reflexLength, reflexLength);
		reflection3D(regCurrentSegmentation, height - reflexLength, length - reflexLength, width - reflexLength, reflexLength);
		//==============================================================================
		// Create a dist. map from the registered bin map
		dataType** regCurrSegDistMap = (dataType**)malloc(sizeof(dataType**) * height);
		for (size_t j = 0; j < height; j++)
		{
			regCurrSegDistMap[j] = (dataType*)malloc(sizeof(dataType*) * (length * width));
		}
		fastSweepingFunction_3D(regCurrSegDistMap, regCurrentSegmentation, length, width, height, 1, 50000, foreground);
		//==============================================================================
		// Signed dist. map
		ReDistFn(regCurrSegDistMap, height, length, width, foreground, background);
		//==============================================================================
		// Default Estimation params
		dataType bound = 3.0, eps = 1.0;
		int est_iterations = 500;
		// Inputs are dist. Maps for Mean shape, current registered segmentation, current segmentation
		if ((*pcaParam).estMethod == MIN_Energy) // Minimizing l2-loss energy
		{
			bound = 1.25, eps = 1.0;
			est_iterations = 200;
			estimatedShapeByEnergy(dta3D,(*pcaParam).dtaMeanDist, regCurrSegDistMap, currentSegmentation, (*pcaParam).eigenvalues, (*pcaParam).eigenvectors, (*pcaParam).princomp, height, length, width, regParams.affineResults, bound, eps, est_iterations, false, 2, h, tolerance);
		}
		else if ((*pcaParam).estMethod == MIN_Probability) // Minimizing probability
		{
			bound = 0.25, eps = 1.0;
			estimateShapeByProbability(dta3D ,(*pcaParam).dtaMeanDist, regCurrSegDistMap, currentSegmentation, (*pcaParam).eigenvalues, (*pcaParam).eigenvectors, (*pcaParam).princomp, height, length, width, regParams.affineResults, bound, eps, est_iterations, false, 2, h, tolerance);
		}
		else
		{
			printf("%s\n", st);
			printf("No estimation method defined, exiting\n");
			printf("%s\n\n", st);
			goto ExitCond;
		}
		//==============================================================================
		// Free Memory - current regitered segmentation
		for (k = 0; k < height + 2 * reflexLength + 1; k++)
		{
			free(regCurrentSegmentation[k]);
		}
		for (size_t j = 0; j < height; j++)
		{
			free(regCurrSegDistMap[j]);
		}
		free(regCurrentSegmentation);
		free(regCurrSegDistMap);
		regCurrSegDistMap = NULL;
		regCurrentSegmentation = NULL;
		//==============================================================================
		reflection3D((*dta3D).priorShape, height - reflexLength, length - reflexLength, width - reflexLength, reflexLength);
		//==============================================================================
		// Smooth Estimated to the same level as current segmentation
		stepsize = 0.05, time_steps = (size_t)smooth_level * 2.5, tolerance = 6.0 * tol_s;
		// Error btn Estimated segmentation and current segmentation
		error = errorCalc((*dta3D).priorShape, currentSegmentation, height, length, width, 1);
		//==============================================================================
		printf("Error between current and estimated segmentation: %.8lf\n", error);
		//==============================================================================
		// MHD
		// Calaulation of hausdorf distance
		double val_a = hausdorffDist(currentSegmentation, (*dta3D).priorShape, height, length, width);
		double val_b = hausdorffDist((*dta3D).priorShape, currentSegmentation, height, length, width);
		double hdistance = max(val_a, val_b);
		printf("Hausdorf distance between current and estimated segmentation: %lf\n", hdistance);
		//==============================================================================
		goto ExitCond;
	}
ExitCond:
	printf("%s\n", st);
	printf("End of estimation of the current segmentation shape\n");
	printf("%s\n\n", st);
	// Return error
	return error;
}
//==============================================================================
void segmentation3D_Ap_coef(ABSContainer* dta3D, GradData* grad3D, AtlasData* atls3D, size_t height, size_t length, size_t width, dataType** priorShape, dataType lambda, dataType eta, dataType zeta, dataType h, size_t p, dataType epsilon, dataType tau, dataType mu1, dataType mu2, dataType d1, dataType d2)
{
	//==============================================================================
	size_t k, i, j, xd;
	dataType uqx, uqy, uqz, Gn, Ge, Gw, Gs, Gt, Gb, h3 = h * h * h, tmp, tmp1, tmp2, tmp3, tmp4, tmp6, tmp7, tmp8, tmp9, tmp10, gama;
	dataType Gpe, Gpw, Gpb, Gpt, Gps, Gpn;
	dataType Apq;
	// Calculate the central data - x-forward,y-backward,z-vertical directions
	for (k = p; k < height; k++)
	{
		for (i = p; i < length; i++)
		{
			for (j = p; j < width; j++)
			{
				// 1D Conversion of row and column
				xd = x_new(i, j, length);
				// Calculations
				(*dta3D).dta_h[k][xd] = ((*dta3D).dta_u[k][x_new(i - 1, j - 1, length)] + (*dta3D).dta_u[k][x_new(i, j - 1, length)] + (*dta3D).dta_u[k][xd] + (*dta3D).dta_u[k][x_new(i - 1, j, length)]) / 4.; // Vertical
				(*dta3D).dta_l[k][xd] = ((*dta3D).dta_u[k - 1][x_new(i - 1, j, length)] + (*dta3D).dta_u[k - 1][xd] + (*dta3D).dta_u[k][xd] + (*dta3D).dta_u[k][x_new(i - 1, j, length)]) / 4.; // Forwardu
				(*dta3D).dta_w[k][xd] = ((*dta3D).dta_u[k - 1][x_new(i, j - 1, length)] + (*dta3D).dta_u[k - 1][xd] + (*dta3D).dta_u[k][xd] + (*dta3D).dta_u[k][x_new(i, j - 1, length)]) / 4.; // Backward
			}
		}
	}
	// Evaluate the gradients for the calculated data
	for (k = p; k < height; k++)
	{
		for (i = p; i < length; i++)
		{
			for (j = p; j < width; j++)
			{
				// 1D Conversion of row and column
				xd = x_new(i, j, length);

				// East
				uqx = ((*dta3D).dta_u[k][x_new(i, j + 1, length)] - (*dta3D).dta_u[k][xd]) / h;
				uqy = ((*dta3D).dta_h[k][x_new(i + 1, j + 1, length)] - (*dta3D).dta_h[k][x_new(i, j + 1, length)]) / h; // y-direction
				uqz = ((*dta3D).dta_w[k + 1][x_new(i, j + 1, length)] - (*dta3D).dta_w[k][x_new(i, j + 1, length)]) / h; // z-direction
				Ge = uqx * uqx + uqy * uqy + uqz * uqz + epsilon;
				// West
				uqx = ((*dta3D).dta_u[k][xd] - (*dta3D).dta_u[k][x_new(i, j - 1, length)]) / h;
				uqy = (((*dta3D).dta_h[k][x_new(i + 1, j, length)] - (*dta3D).dta_h[k][xd]) / h); // y-direction
				uqz = (((*dta3D).dta_w[k + 1][xd] - (*dta3D).dta_w[k][xd]) / h); // z-direction
				Gw = uqx * uqx + uqy * uqy + uqz * uqz + epsilon;
				// North
				uqy = ((*dta3D).dta_u[k][xd] - (*dta3D).dta_u[k][x_new(i - 1, j, length)]) / h; // y-direction
				uqx = ((*dta3D).dta_h[k][x_new(i, j + 1, length)] - (*dta3D).dta_h[k][xd]) / h; // x-direction
				uqz = ((*dta3D).dta_l[k + 1][xd] - (*dta3D).dta_l[k][xd]) / h; // z-direction
				Gn = uqx * uqx + uqy * uqy + uqz * uqz + epsilon;
				// South
				uqy = (((*dta3D).dta_u[k][x_new(i + 1, j, length)] - (*dta3D).dta_u[k][xd]) / h); // y-direction
				uqx = ((*dta3D).dta_h[k][x_new(i + 1, j + 1, length)] - (*dta3D).dta_h[k][x_new(i + 1, j, length)]) / h; // x-direction
				uqz = ((*dta3D).dta_l[k + 1][x_new(i + 1, j, length)] - (*dta3D).dta_l[k][x_new(i + 1, j, length)]) / h; //z-direction
				Gs = uqx * uqx + uqy * uqy + uqz * uqz + epsilon;
				// Bottom
				uqz = ((*dta3D).dta_u[k][xd] - (*dta3D).dta_u[k - 1][xd]) / h; // z-direction
				uqy = ((*dta3D).dta_l[k][x_new(i + 1, j, length)] - (*dta3D).dta_l[k][xd]) / h; // y-direction
				uqx = ((*dta3D).dta_w[k][x_new(i, j + 1, length)] - (*dta3D).dta_w[k][xd]) / h; // x-direction
				Gb = uqx * uqx + uqy * uqy + uqz * uqz + epsilon;
				// Top
				uqz = ((*dta3D).dta_u[k + 1][xd] - (*dta3D).dta_u[k][xd]) / h; // z-direction
				uqy = ((*dta3D).dta_l[k + 1][x_new(i + 1, j, length)] - (*dta3D).dta_l[k + 1][xd]) / h; // y-direction
				uqx = ((*dta3D).dta_w[k + 1][x_new(i, j + 1, length)] - (*dta3D).dta_w[k + 1][xd]) / h; // x-direction
				Gt = uqx * uqx + uqy * uqy + uqz * uqz + epsilon;
				//==============================================================================
				// Approximate gradients of edge detector function
				Gpe = ((*grad3D).gd[k][x_new(i, j + 1, length)] - (*grad3D).gd[k][x_new(i, j - 1, length)]) / (2 * h);
				Gpw = -1 * Gpe;
				Gpn = ((*grad3D).gd[k][x_new(i - 1, j, length)] - (*grad3D).gd[k][x_new(i + 1, j, length)]) / (2 * h);
				Gps = -1 * Gpn;
				Gpt = ((*grad3D).gd[k + 1][xd] - (*grad3D).gd[k - 1][xd]) / (2 * h);
				Gpb = -1 * Gpt;
				//==============================================================================
				/*
				Gpe = (*grad3D).ge[k][xd] - (*grad3D).gw[k][xd];
				Gpw = -1 * Gpe;
				Gpn = (*grad3D).gn[k][xd] - (*grad3D).gs[k][xd];
				Gps = -1 * Gpn;
				Gpt = (*grad3D).gt[k][xd] - (*grad3D).gb[k][xd];
				Gpb = -1 * Gpt;
				*/
				//==============================================================================
				//if (pcaParam == NULL)
				//{
				//	printf("Exiting... NULL pca parameters\n");
				//	return;
				//}
				//==============================================================================
				//tmp4 = mu3 * (*dta3D).dta_u[k][xd] - priorShape[k][xd]; // u - u`
				//gama = GammaFun(tmp4); // |u - u|
				//if (gama != gama)
				//{
				//	printf("nan gama\n");
				//}
				//printf("\n The absolute diffenceis:%.f\n", fabsf(tmp4));
				//printf("\n The value of gamma is:%.f\n", gama);
				//==============================================================================
				// tmp holder
				tmp = (tau / h3) * sqrt((Gw + Ge + Gs + Gn + Gb + Gt) / 6.0);
				if (tmp != tmp)
				{
					printf("nan tmp\n");
				}
				// Calculated data gradients
				(*grad3D).aw[k][xd] = tmp * g_mcf(Gw) * (*grad3D).gw[k][xd] * mu2;
				(*grad3D).ae[k][xd] = tmp * g_mcf(Ge) * (*grad3D).ge[k][xd] * mu2;
				(*grad3D).as[k][xd] = tmp * g_mcf(Gs) * (*grad3D).gs[k][xd] * mu2;
				(*grad3D).an[k][xd] = tmp * g_mcf(Gn) * (*grad3D).gn[k][xd] * mu2;
				(*grad3D).ab[k][xd] = tmp * g_mcf(Gb) * (*grad3D).gb[k][xd] * mu2;
				(*grad3D).at[k][xd] = tmp * g_mcf(Gt) * (*grad3D).gt[k][xd] * mu2;
				tmp1 = ((*grad3D).aw[k][xd] + (*grad3D).ae[k][xd] + (*grad3D).as[k][xd] + (*grad3D).an[k][xd] + (*grad3D).at[k][xd] + (*grad3D).ab[k][xd]);
				if (tmp1 != tmp1)
				{
					printf("nan tmp1\n");
				}
				//==============================================================================
				tmp2 = (tau / h3);
				//==============================================================================

				//tmp6 = (((*dta3D).dta_u[k][xd] - priorShape[k][xd]) + ((*dta3D).dta_u[k][x_new(i, j - 1, length)] - priorShape[k][x_new(i, j - 1, length)])) / 2.;
				tmp6 = (((*dta3D).dta_u[k][xd] + (*dta3D).dta_u[k][x_new(i, j - 1, length)]) - (priorShape[k][xd] + priorShape[k][x_new(i, j - 1, length)])) / 2.;
				//tmp6 = ((*dta3D).dta_u[k][x_new(i, j - 1, length)] - priorShape[k][x_new(i, j - 1, length)]);
				gama = GammaFun(tmp6, d1, d2); // |u - u|

				tmp9 = ((*dta3D).dta_u[k][x_new(i, j - 1, length)] - (*dta3D).dta_u[k][xd]) * (g_mcf(Gw));
				tmp7 = gama * (1 - lambda) * (zeta) * ((((*atls3D).gd[k][x_new(i, j - 1, length)] + (*atls3D).gd[k][xd]) / 2) * (tmp9));
				tmp8 = (gama * lambda * Gpw) * (eta);

				tmp10 = (1 - gama) * (tmp6 * tmp9);

				Apq = tmp7 - tmp8 + tmp10;
				(*atls3D).bp_w[k][xd] = mu1 * tmp2 * min(Apq, 0);
				//==============================================================================
				//tmp6 = (((*dta3D).dta_u[k][xd] - priorShape[k][xd]) + ((*dta3D).dta_u[k][x_new(i - 1, j, length)] - priorShape[k][x_new(i - 1, j, length)])) / 2;
				tmp6 = (((*dta3D).dta_u[k][xd] + (*dta3D).dta_u[k][x_new(i - 1, j, length)]) - (priorShape[k][xd] + priorShape[k][x_new(i - 1, j, length)])) / 2;
				//tmp6 = ((*dta3D).dta_u[k][x_new(i - 1, j, length)] - priorShape[k][x_new(i - 1, j, length)]);
				gama = GammaFun(tmp6, d1, d2); // |u - u|

				tmp9 = ((*dta3D).dta_u[k][x_new(i - 1, j, length)] - (*dta3D).dta_u[k][xd]) * (g_mcf(Gn));
				tmp7 = gama * (1 - lambda) * (zeta) * ((((*atls3D).gd[k][x_new(i - 1, j, length)] + (*atls3D).gd[k][xd]) / 2) * (tmp9));
				tmp8 = (gama * lambda * Gpn) * (eta);

				tmp10 = (1 - gama) * (tmp6 * tmp9);

				Apq = tmp7 - tmp8 + tmp10;
				(*atls3D).bp_n[k][xd] = mu1 * tmp2 * min(Apq, 0);
				//==============================================================================
				//tmp6 = (((*dta3D).dta_u[k][xd] - priorShape[k][xd]) + ((*dta3D).dta_u[k][x_new(i + 1, j, length)] - priorShape[k][x_new(i + 1, j, length)])) / 2;
				tmp6 = (((*dta3D).dta_u[k][xd] + (*dta3D).dta_u[k][x_new(i + 1, j, length)]) - (priorShape[k][xd] + priorShape[k][x_new(i + 1, j, length)])) / 2;
				//tmp6 = ((*dta3D).dta_u[k][x_new(i + 1, j, length)] - priorShape[k][x_new(i + 1, j, length)]);
				gama = GammaFun(tmp6, d1, d2); // |u - u|

				tmp9 = ((*dta3D).dta_u[k][x_new(i + 1, j, length)] - (*dta3D).dta_u[k][xd]) * (g_mcf(Gs));
				tmp7 = gama * (1 - lambda) * (zeta) * ((((*atls3D).gd[k][x_new(i + 1, j, length)] + (*atls3D).gd[k][xd]) / 2) * (tmp9));
				tmp8 = (gama * lambda * Gps) * (eta);

				tmp10 = (1 - gama) * (tmp6 * tmp9);

				Apq = tmp7 - tmp8 + tmp10;
				(*atls3D).bp_s[k][xd] = mu1 * tmp2 * min(Apq, 0);
				//==============================================================================
				//tmp6 = (((*dta3D).dta_u[k][xd] - priorShape[k][xd]) + ((*dta3D).dta_u[k + 1][xd] - priorShape[k + 1][xd])) / 2;
				tmp6 = (((*dta3D).dta_u[k][xd] + (*dta3D).dta_u[k + 1][xd]) - (priorShape[k][xd] + priorShape[k + 1][xd])) / 2;
				//tmp6 = ((*dta3D).dta_u[k + 1][xd] - priorShape[k + 1][xd]);
				gama = GammaFun(tmp6, d1, d2); // |u - u|

				tmp9 = ((*dta3D).dta_u[k + 1][xd] - (*dta3D).dta_u[k][xd]) * g_mcf(Gt);
				tmp7 = gama * (1 - lambda) * (zeta) * ((((*atls3D).gd[k + 1][xd] + (*atls3D).gd[k][xd]) / 2) * (tmp9));
				tmp8 = (gama * lambda * Gpt) * (eta);

				tmp10 = (1 - gama) * (tmp6 * tmp9);

				Apq = tmp7 - tmp8 + tmp10;
				(*atls3D).bp_t[k][xd] = mu1 * tmp2 * min(Apq, 0);
				//==============================================================================
				//tmp6 = (((*dta3D).dta_u[k][xd] - priorShape[k][xd]) + ((*dta3D).dta_u[k - 1][xd] - priorShape[k - 1][xd])) / 2;
				tmp6 = (((*dta3D).dta_u[k][xd] + (*dta3D).dta_u[k - 1][xd]) - (priorShape[k][xd] + priorShape[k - 1][xd])) / 2;
				//tmp6 = ((*dta3D).dta_u[k - 1][xd] - priorShape[k - 1][xd]);
				gama = GammaFun(tmp6, d1, d2); // |u - u|

				tmp9 = ((*dta3D).dta_u[k - 1][xd] - (*dta3D).dta_u[k][xd]) * g_mcf(Gb);
				tmp7 = gama * (1 - lambda) * (zeta) * ((((*atls3D).gd[k - 1][xd] + (*atls3D).gd[k][xd]) / 2) * (tmp9));
				tmp8 = (gama * lambda * Gpb) * (eta);

				tmp10 = (1 - gama) * (tmp6 * tmp9);

				Apq = tmp7 - tmp8 + tmp10;
				(*atls3D).bp_b[k][xd] = mu1 * tmp2 * min(Apq, 0);
				//==============================================================================
				//tmp6 = (((*dta3D).dta_u[k][xd] - priorShape[k][xd]) + ((*dta3D).dta_u[k][x_new(i, j + 1, length)] - priorShape[k][x_new(i, j + 1, length)])) / 2;
				tmp6 = (((*dta3D).dta_u[k][xd] + (*dta3D).dta_u[k][x_new(i, j + 1, length)]) - (priorShape[k][xd] + priorShape[k][x_new(i, j + 1, length)])) / 2;
				//tmp6 = ((*dta3D).dta_u[k][x_new(i, j + 1, length)] - priorShape[k][x_new(i, j + 1, length)]);
				gama = GammaFun(tmp6, d1, d2); // |u - u|

				tmp9 = ((*dta3D).dta_u[k][x_new(i, j + 1, length)] - (*dta3D).dta_u[k][xd]) * g_mcf(Ge);
				tmp7 = gama * (1 - lambda) * (zeta) * ((((*atls3D).gd[k][x_new(i, j + 1, length)] + (*atls3D).gd[k][xd]) / 2) * (tmp9));
				tmp8 = (gama * lambda * Gpe) * (eta);

				tmp10 = (1 - gama) * (tmp6 * tmp9);

				Apq = tmp7 - tmp8 + tmp10;
				(*atls3D).bp_e[k][xd] = mu1 * tmp2 * min(Apq, 0);
				//==============================================================================
				//tmp3 = ((*atls3D).bp_e[k][xd] + (*atls3D).bp_w[k][xd] + (*atls3D).bp_n[k][xd] + (*atls3D).bp_s[k][xd] + (*atls3D).bp_t[k][xd] + (*atls3D).bp_b[k][xd]);
				//==============================================================================
				// ap
				//(*grad3D).ap[k][xd] = 1 - tmp3 + tmp1;
				(*grad3D).ap[k][xd] = 1 + tmp1;
				if ((*grad3D).ap[k][xd] == 0)
				{
					printf("(*grad3D).ap[k][xd] equals 0\n");
				}
				//==============================================================================
			}
		}
	}
}
//==============================================================================
void gmcf3D_atlas(ABSContainer* dta3D, GradData* grad3D, AtlasData* atls3D, size_t height, size_t length, size_t width, size_t p)
{
	size_t k, i, j, xd;
	int z, ni = 2000, count;
	dataType y, error, tmp, tmp2, tmp3;;
	//==============================================================================
	dataType omega_s = 1.12, tol_s = 1.0e-06;
	size_t smooth_level = 750;
	//==============================================================================
	dataType imgFg = 0.0f, imgBg = 1.0f;
	// Reset U0 values before copying in to
	resetDta((*dta3D).dta_u0, height, length, width, imgFg);
	//==============================================================================
	// Copy u to u0
	for (k = 0; k < height; k++)
	{
		for (i = 0; i < length; i++)
		{
			for (j = 0; j < width; j++)
			{
				// 1D Conversion of row and column
				xd = x_new(i, j, length);

				(*dta3D).dta_u0[k][xd] = (*dta3D).dta_u[k][xd];
			}
		}
	}
	z = 0;
	do
	{
		z = z + 1;
		for (k = 0; k < height; k++)
		{
			for (i = 0; i < length; i++)
			{
				for (j = 0; j < width; j++)
				{
					// 1D Conversion of row and column
					xd = x_new(i, j, length);
					// Begin Gauss Seidel SOR
					//==============================================================================
					tmp2 = ((*grad3D).aw[k][xd] * (*dta3D).dta_u[k][x_new(i, j - 1, length)]
						+ (*grad3D).ae[k][xd] * (*dta3D).dta_u[k][x_new(i, j + 1, length)]
						+ (*grad3D).as[k][xd] * (*dta3D).dta_u[k][x_new(i + 1, j, length)]
						+ (*grad3D).an[k][xd] * (*dta3D).dta_u[k][x_new(i - 1, j, length)]
						+ (*grad3D).at[k][xd] * (*dta3D).dta_u[k + 1][xd]
						+ (*grad3D).ab[k][xd] * (*dta3D).dta_u[k - 1][xd]);
					if (tmp2 != tmp2)
					{
						printf("nan tmp2\n");
					}
					//==============================================================================
					tmp3 = ((*atls3D).bp_w[k][xd] * ((*dta3D).dta_u[k][x_new(i, j - 1, length)] - (*dta3D).dta_u[k][xd])
						+ (*atls3D).bp_e[k][xd] * ((*dta3D).dta_u[k][x_new(i, j + 1, length)] - (*dta3D).dta_u[k][xd])
						+ (*atls3D).bp_s[k][xd] * ((*dta3D).dta_u[k][x_new(i + 1, j, length)] - (*dta3D).dta_u[k][xd])
						+ (*atls3D).bp_n[k][xd] * ((*dta3D).dta_u[k][x_new(i - 1, j, length)] - (*dta3D).dta_u[k][xd])
						+ (*atls3D).bp_t[k][xd] * ((*dta3D).dta_u[k + 1][xd] - (*dta3D).dta_u[k][xd])
						+ (*atls3D).bp_b[k][xd] * ((*dta3D).dta_u[k - 1][xd] - (*dta3D).dta_u[k][xd]));
					if (tmp3 != tmp3)
					{
						printf("nan tmp3\n");
					}
					//==============================================================================
					y = ((*dta3D).dta_u0[k][xd] - tmp3 + tmp2) / (*grad3D).ap[k][xd];
					if (isinf(y))
					{
						printf("inf y\n");
					}
					// Apply the SOR y
					if (y != y)
					{
						printf("nan y\n");
					}
					(*dta3D).dta_u[k][xd] = (*dta3D).dta_u[k][xd] + omega_s * (y - (*dta3D).dta_u[k][xd]);
				}
			}
		}
		//==============================================================================
		//resetDta((*dta3D).dta_tmp, foregound);
		//fastSweepingFunction_3D((*dta3D).dta_tmp, (*dta3D).dta_u, length, width, height, 1, 50000, foregound);
		//==============================================================================
		error = 0.0;
		count = 0;
		for (k = p; k < height; k++)
		{
			for (i = p; i < length; i++)
			{
				for (j = p; j < width; j++)
				{
					// 1D Conversion of row and column
					xd = x_new(i, j, length);
					// Error calculation
					tmp2 = ((*grad3D).aw[k][xd] * (*dta3D).dta_u[k][x_new(i, j - 1, length)]
						+ (*grad3D).ae[k][xd] * (*dta3D).dta_u[k][x_new(i, j + 1, length)]
						+ (*grad3D).as[k][xd] * (*dta3D).dta_u[k][x_new(i + 1, j, length)]
						+ (*grad3D).an[k][xd] * (*dta3D).dta_u[k][x_new(i - 1, j, length)]
						+ (*grad3D).at[k][xd] * (*dta3D).dta_u[k + 1][xd]
						+ (*grad3D).ab[k][xd] * (*dta3D).dta_u[k - 1][xd]);
					if (tmp2 != tmp2)
					{
						printf("nan tmp2\n");
					}
					//==============================================================================
					tmp3 = ((*atls3D).bp_w[k][xd] * (*dta3D).dta_u[k][x_new(i, j - 1, length)]
						+ (*atls3D).bp_e[k][xd] * (*dta3D).dta_u[k][x_new(i, j + 1, length)]
						+ (*atls3D).bp_s[k][xd] * (*dta3D).dta_u[k][x_new(i + 1, j, length)]
						+ (*atls3D).bp_n[k][xd] * (*dta3D).dta_u[k][x_new(i - 1, j, length)]
						+ (*atls3D).bp_t[k][xd] * (*dta3D).dta_u[k + 1][xd]
						+ (*atls3D).bp_b[k][xd] * (*dta3D).dta_u[k - 1][xd]);
					if (tmp3 != tmp3)
					{
						printf("nan tmp3\n");
					}
					//==============================================================================
					tmp = (*grad3D).ap[k][xd] * (*dta3D).dta_u[k][xd] - (*dta3D).dta_u0[k][xd] + tmp3 - tmp2;
					if (tmp != tmp)
					{
						printf("nan tmp\n");
					}
					error += (tmp * tmp);
					if (error != error)
					{
						printf("nan error\n");
					}
					count++;
				}
			}
		}
		// Mean Square Error
		if (count == 0)
		{
			count = 1;
		}
		error = error / count;
	} while ((error > tol_s) && (z < ni));
	printf("The number of GMCF gaus-seidel calulation iterations is %d\n", z);
	printf("Error is %e\n", error);
	//==============================================================================
	// Redistance the seg fn
	ReDistFn((*dta3D).dta_u, height, length, width, imgFg, imgBg);
	//==============================================================================
	dataType stepsize = 0.05, time_steps = smooth_level, tolerance = 6.0 * tol_s;
	//smoothData((*dta3D).dta_u, height, length, width, stepsize, time_steps, tolerance);
	//==============================================================================
}
//==============================================================================
dataType calcRho(dataType** aPtr, dataType** bPtr, dataType distWithin, size_t height, size_t length, size_t width, size_t p) {
	size_t k, i, j, xd;
	dataType rho = 0;
	int counts = 0;
	for (k = p; k < height; k++)
	{
		for (i = p; i < length; i++)
		{
			for (j = p; j < width; j++)
			{
				// 1D Conversion of row and column
				xd = x_new(i, j, length);

				if (bPtr[k][xd] < distWithin) // Inside the region defined for init seg fn
				{
					rho += aPtr[k][xd]; // adds intensity in original image to rho
					counts++;
				}
			}
		}
	}
	// Find average intesity
	if (counts == 0) // Ensuring no division by zero!
	{
		counts = 1;
	}
	rho = rho / counts;
	return rho;
}
//==============================================================================
void calcVariance(dataType** aPtr, dataType** bPtr, dataType distWithin, size_t height, size_t length, size_t width, size_t p, int step, dataType rho_zero_step) {
	size_t k, i, j, xd;
	// Calc. rho
	dataType rho = calcRho(aPtr, bPtr, distWithin, height, length, width, p);
	dataType variance = 0.;
	int counts = 0;
	// Calc the sample variance
	for (k = p; k < height; k++)
	{
		for (i = p; i < length; i++)
		{
			for (j = p; j < width; j++)
			{
				// 1D Conversion of row and column
				xd = x_new(i, j, length);

				if (bPtr[k][xd] < distWithin) // Inside the region
				{
					//T tmp = aPtr[k][xd] - rho;
					dataType tmp = aPtr[k][xd] - rho_zero_step;
					variance += tmp * tmp;
					counts++;
				}
			}
		}
	}
	if (counts == 0) // Ensuring no division by zero!
	{
		counts = 1;
	}
	variance /= (counts - 1);
	//==============================================================================
}
//==============================================================================
dataType error_3D(dataType** aPtr, dataType** bPtr, size_t height, size_t length, size_t width, dataType h_val, size_t step) {
	size_t i, j, k, xd;
	dataType** tmPtr = (dataType**)malloc(sizeof(dataType*) * height);
	const size_t dim2D_mem = sizeof(dataType) * (length * width);
	for (i = 0; i < height; i++)
	{
		tmPtr[i] = (dataType*)malloc(dim2D_mem);
	}
	dataType tmp, tmax = 1.0e+38 * -1.0;
	for (k = 0; k < height; k++)
	{
		for (i = 0; i < length; i++)
		{
			for (j = 0; j < width; j++)
			{
				// 1D Conversion of row and column
				xd = x_new(i, j, length);
				// Error calculation
				tmp = (aPtr[k][xd] - bPtr[k][xd]) * h_val;
				tmp = fabs(tmp);
				// Fill ptr
				tmPtr[k][xd] = tmp;
				// find max diff.
				if (tmp > tmax)
				{
					tmax = tmp;
				}
			}
		}
	}
	//==============================================================================
	for (i = 0; i < height; i++)
	{
		free(tmPtr[i]);
	}
	free(tmPtr);
	tmPtr = NULL;
	//==============================================================================
	return tmax;
}
//==============================================================================
// Curve estimation functions

// Extend array function
void extend_array(dataType** ext_array, dataType* in_array, int length, int w_size)
{
	int i;
	*ext_array = (dataType*)malloc(sizeof(dataType) * (length + w_size));;
	// Copy to extended array
	int count_ext = 2;
	for (i = 0; i < (length + w_size); i++)
	{
		if (i < length)
		{
			(*ext_array)[i] = in_array[i];
		}
		else if (i >= length)
		{
			(*ext_array)[i] = in_array[length - count_ext];
			count_ext++;
		}
	}
}
// Moving average fn
void moving_average(dataType** results, dataType* dtarray, int w_size, int length)
{
	int i;
	*results = (dataType*)malloc(sizeof(dataType) * (length));
	// Extend
	dataType* ext_arr = NULL;
	extend_array(&ext_arr, dtarray, length, w_size);
	dataType tmp;
	for (i = 0; i < length; i++)
	{
		tmp = (ext_arr[i] + ext_arr[i + 1] + ext_arr[i + 2] + ext_arr[i + 3] + ext_arr[i + 4]) / (dataType)w_size;
		(*results)[i] = tmp;
	}
	// Free
	free(ext_arr);
	ext_arr = NULL;
}
// Est. a curve by offset val.
void est_curve(dataType** est_y, dataType* arr, dataType* d_arr, int count_b, int count_e, int* est_x, int length, int offset)
{
	int i;
	// array for idx + 1 of array
	dataType* orig_y = (dataType*)malloc(sizeof(dataType) * (length - 1));
	int* orig_x = malloc(sizeof(int) * (length - 1));
	*est_y = (dataType*)malloc(sizeof(dataType) * (length - 1));
	for (i = 0; i < length - 1; i++)
	{
		(*est_y)[i] = 0.;
		est_x[i] = 0;
	}
	for (i = 0; i < length; i++)
	{
		orig_y[i] = arr[i + 1];
		orig_x[i] = i + 1;
	}
	for (i = count_b; i < count_e; i++)
	{
		dataType dy = -0.2 * (d_arr[i - 2] + d_arr[i - 1] + d_arr[i] + d_arr[i + 1] + d_arr[i + 2]);
		est_x[i] = orig_x[i] + offset;
		(*est_y)[i] = orig_y[i] + (dy * offset);
		if ((*est_y)[i] < 0)
		{
			(*est_y)[i] = 0.;
		}
	}
	// Free
	/*free(orig_y);
	free(orig_x);*/
	orig_x = NULL;
	orig_y = NULL;
}
// Area of curve
dataType curve_area(dataType* arr, int length, int interval)
{
	int i;
	dataType total = 0.;
	for (i = 0; i < length; i++)
	{
		if ((i == 0) || (i == length - 1))
		{
			total += arr[i] / 2.;
		}
		else
		{
			total += arr[i];
		}
	}
	return total * interval;
}
// Half area for a curve
int half_area(dataType* arr, int count_a, int length)
{
	int i;
	dataType* f_area = (dataType*)malloc(sizeof(dataType) * (length));
	// Initialize to zero
	for (i = 0; i < length; i++)
	{
		// Initial area
		if (i < 2)
		{
			f_area[i] = arr[i];
		}
		else
		{
			f_area[i] = 0.;
		}
	}
	// Initial area
	dataType h_area = 0.5 * curve_area(arr, length, 1);
	for (i = count_a; i < length + 1; i++)
	{
		dataType _area = curve_area(arr, i, 1);
		if (_area >= h_area)
		{
			printf("Half area %.6lf at index %d\n", _area, i);
			break;
		}
		f_area[i] = arr[i];
	}
	// Free
	free(f_area);
	return i;
}
// Find min
int find_min(dataType* arr, dataType* minimum, int begin, int length)
{
	int index, count = 0;
	(*minimum) = arr[begin];
	int i;
	for (i = begin; i < length; i++)
	{
		if ((arr[i] < (*minimum)) && (i != begin))
		{
			(*minimum) = arr[i];
			index = i;
		}
		if (((*minimum) == 0) && (i != begin))
		{
			count = count + begin;
			break;
		}
		count++;
	}
	return i;
}
// Difference btn idx and idx + 1 array values
void diff_s(dataType** rlts, dataType* arr, int length)
{
	int i;
	*rlts = (dataType*)malloc(sizeof(dataType) * (length - 1));
	for (i = 0; i < length - 1; i++)
	{
		(*rlts)[i] = arr[i] - arr[i + 1];
	}
}
//==============================================================================
// Stop segmentation function
bool stop_segment3D(
	ABSContainer* dta3D,
	tmpDataHolders* tmpDataStepHolder,
	size_t p, dataType h3, size_t* step,
	double* lambda, double* zeta, double* eta,
	dataType* mass_prev, size_t* max_iters,
	dataType tol_m, dataType tol_e,
	PCAData* pcaParam,Optimization_Method optMethod, Registration_Params regParams
) {
	//size_t i, j;
	dataType mass = 0.0;
	bool stpcond = false;
	//==============================================================================
	dataType foreground = (*dta3D).imgForeground, background = (*dta3D).imgBackground;
	//==============================================================================
	size_t height = (*dta3D).height, length = (*dta3D).length, width = (*dta3D).width;
	//==============================================================================
	// calc rho, variance
	//calcVar((*dta3D).dta, (*dta3D).dta_u, 0.7, height, length, width, p, P, (*step), statsPtr, tmpDataStepHolder->rho_zero_step);
	calcVar((*dta3D).dta, (*dta3D).dta_u, 0.0, height, length, width, p, (*step), tmpDataStepHolder->rho_zero_step);
	//==============================================================================
	// Apply minimization of least square
	// https://en.wikipedia.org/wiki/Stochastic_gradient_descent
	//==============================================================================
	mass = errorCalc((*dta3D).dta_u0, (*dta3D).dta_u, height, length, width, h3);
	//==============================================================================
	double diff_mass = fabs(mass - (*mass_prev));
	//==============================================================================
	// Update previous weight after
	(*mass_prev) = mass;
	//==============================================================================
	// Save the 3D difference to a file
	dataType max_diff = error_3D((*dta3D).dta_u0, (*dta3D).dta_u, height, length, width, h3, (*step));
	//==============================================================================
	// Useful pointers for below fn's
	int begin = 2, end;
	//==============================================================================
	// Determine what step to go back and reduce g2 effect and also where to stop it
	// 0. Store the mass in a 1D array
	// 1. smooth fn
	// 2. Difference
	// 3. Estimate curve and find first zero
	// 4. If 1st zero, determin how far to go back - 10th or 20th step
	// 5 Stop the above
	if (tmpDataStepHolder->est_fun) // Only do this when still estimating
	{
		// Save the mass vals
		tmpDataStepHolder->d_mass[(*step) - 1] = mass;
		if ((*step) == 2)
		{
			tmpDataStepHolder->d_mass[(*step) - 2] = tmpDataStepHolder->d_mass[(*step) - 1];
			//tmpDataStepHolder->d_mass[(*step) - 1] = 0.;
		}
		// Start smoothing the results after steps > 10
		if ((*step) > tmpDataStepHolder->begie_est)
		{
			// Smooth the current results
			dataType* smooth_ptrs = NULL;
			moving_average(&smooth_ptrs, tmpDataStepHolder->d_mass, tmpDataStepHolder->w_size, (*step));
			// The differences
			dataType* diff = NULL;
			diff_s(&diff, smooth_ptrs, (*step));
			// Estimate the current curve
			int* est_x = (int*)malloc(sizeof(int) * ((*step) - 1));
			end = (*step) - begin - 1;
			dataType* est_y = NULL;
			est_curve(&est_y, smooth_ptrs, diff, begin, end, est_x, (*step), tmpDataStepHolder->offset);
			// Finding minimum
			dataType minimum = -1.;
			int idx_min = find_min(est_y, &minimum, begin, end);
			printf("Minimum %lf found at index %d\n", minimum, est_x[idx_min]);
			if ((minimum == 0))
			{
				// Calculate half area, also first zero found
				int _b = idx_min - tmpDataStepHolder->offset + 1, _e = idx_min + 1;
				// Previous 10 values up to the 1st zero.
				// Replace the last 10 values with values from estimated
				for (size_t i = _b; i < _e; i++)
				{
					tmpDataStepHolder->d_mass[i + 1] = est_y[i];
				}
				// Half area
				int idx_hlf_area = half_area(tmpDataStepHolder->d_mass, begin, est_x[idx_min]);
				// Define params to determine skipping g2
				printf("Found first zero at %d\n", est_x[idx_min]);
				// Set where to start reduce g2
				tmpDataStepHolder->beg_reduce = idx_hlf_area;
				//==============================================================================
				// ADjust max iteration with the changes also
				//(*max_iters) = (*max_iters) - (*step);
				//(*max_iters) = (*max_iters) - est_x[idx_min];
				//(*max_iters) = 82;
				//==============================================================================
				// Reset the steps and data to initial
				(*step) = 0;
				// Set to go back to begin reduce step results
				copyDataPointer(tmpDataStepHolder->d_0th_step, (*dta3D).dta_u, height + 2 * p + 1, length + 2 * p + 1, width + 2 * p + 1);
				//==============================================================================
				// Set to false to stop with the estimations
				tmpDataStepHolder->est_fun = false;
				// Activate reducing g2
				tmpDataStepHolder->reduce_g2 = true;
				// Free mem at the end
				free(diff);
				free(smooth_ptrs);
				free(est_x);
				free(est_y);
				smooth_ptrs = NULL;
				diff = NULL;
				est_x = NULL;
				est_y = NULL;
			}
		}
	}
	//==============================================================================
	// G2 skip
	if (!tmpDataStepHolder->est_fun)
	{
		//==============================================================================
		// Only doit once!
		if ((!(*dta3D).est_start))
		{
			//==============================================================================
			//if (((*step) == tmpDataStepHolder->beg_reduce) && !(tmpDataStepHolder->skip_g2))
			if ((*step) == tmpDataStepHolder->beg_reduce)
			{
				printf("Decreasing g2 influence from %d step\n\n", tmpDataStepHolder->beg_reduce);
				//==============================================================================
				// New approach - affects only g2
				(*zeta) = tmpDataStepHolder->reduce_zeta_value * (*zeta);
				(*eta) = tmpDataStepHolder->reduce_eta_value * (*eta);
				printf("New reduced Zeta %.4lf\n\n", (*zeta));
				printf("New reduced eta %.4lf\n\n", (*eta));
				//==============================================================================
				tmpDataStepHolder->skip_g2 = true;
				//==============================================================================
				// Set turn off g2
				tmpDataStepHolder->turnoff_g2 = tmpDataStepHolder->beg_reduce + (*step);
				//==============================================================================
			}
			if (((*step) == tmpDataStepHolder->turnoff_g2) && (tmpDataStepHolder->skip_g2))
			{
				printf("Skipping g2 from step %d\n\n", (*step) + 1);
				//==============================================================================
				// Old approach - affects both g1, g2
				// Set lambda to 1.0 to skip using g2 and use only g1
				(*lambda) = 1.0;
				printf("New lambda %.4lf\n\n", (*lambda));
				//==============================================================================
				// Activate atlas
				(*dta3D).est_start = true;
				//==============================================================================
			}
		}
	}
	//==============================================================================
	if ((mass > tol_m) && ((*step) <= (*max_iters)) && (diff_mass != 0.0))
	{
		//if (((*dta3D).est_start) && (step % begin_est_step == 0)) // After N Steps
		if (((*dta3D).est_start))
		{
			printf("\n Estimating current segmentation at Time Step: %d-th\n", (*step));
			bool stp = false;
			// Make a copy of current segmentation to tmp
			//copyDataPointer((*dta3D).dta_u, (*dta3D).dta_tmp, (height + 2 * p + 1), (length + 2 * p + 1), (width + 2 * p + 1));
			copyDataPointer((*dta3D).dta_u, (*dta3D).dta_tmp, (height), (length), (width));
			// Evaluate estimated segmentation shape error
			if (estmateSegmentation(
				dta3D, pcaParam, (*dta3D).dta_tmp,
				height, length, width,
				p, foreground, background,
				tol_m, h3,
				regParams, optMethod
			) < tol_e)
			{
				printf("\nEstimated shape found very similar to current segmentation at step: %d\n", (*step));
				stpcond = true;
			}
			else
			{
				if ((*step) == (*max_iters))
				{
					goto ExitFn;
				}
				else
				{
					goto Cont;
				}
			}
		}
		else
		{
			if ((*step) == (*max_iters))
			{
				goto ExitFn;
			}
			else
			{
				goto Cont;
			}
		}
	Cont:
		printf(" Calculated Mass: %.12lf\n", mass);
		if (stpcond)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	else
	{
		goto ExitFn;
	}
ExitFn:
	printf("Stop Mass: %.12lf\n", mass);
	stpcond = true;
	goto Cont;
	//==============================================================================
}
//==============================================================================
// Atlas Segmentation Model Interface
void atlasSegmentationModel(Image_Data imageData, Segmentation_Paramereters segParameters, PCAData* pcaParam, tmpDataHolders* tmpDataStepHolder)
{
	//==============================================================================
	// Short varaible names abbreviations
	size_t p = segParameters.reflectionLength;
	dataType mu1 = segParameters.mu1, mu2 = segParameters.mu2; // values used in thesis -> mu1 5.0e-01, mu2 5.0e-01
	dataType d1 = segParameters.d1, d2 = segParameters.d2; // values used in thesis -> d1 20.0 (12.0 for 75% noise, 3.0 for 75% noise)
	double lambda = segParameters.lambda, eta = segParameters.eta, zeta = segParameters.zeta;
	size_t imageHeight = imageData.height, imageLength = imageData.length, imageWidth = imageData.width, dim2D = imageLength * imageLength;
	size_t max_iters = segParameters.timeSteps;
	dataType tol_s = segParameters.toleranceSegmentation;
	dataType tol_e = segParameters.toleranceEstimation;
	dataType h3 = segParameters.hSpacing.x * segParameters.hSpacing.y * segParameters.hSpacing.z;
	dataType epsilon = segParameters.episonRegularization, tau = segParameters.tau;
	Optimization_Method optMethod = segParameters.optMethod;
	Registration_Params regParams = segParameters.regParams;
	//==============================================================================
	// Holder for the operations todo after each time step
	const size_t mem_alloc_h = sizeof(dataType*) * (imageHeight + 2 * p + 1);
	const size_t mem_alloc_lw = ((imageLength + 2 * p + 1) * (imageWidth + 2 * p + 1)) * sizeof(dataType);
	const size_t dimLW = ((imageLength + 2 * p + 1) * (imageWidth + 2 * p + 1));
	//==============================================================================
	// Set up data containers for the Gradient data and Atlas Data
	// Atlas
	AtlasData* atls3D;
	initializeAtlasData(atls3D, imageHeight, imageLength, imageWidth, p);
	// Gradient
	GradData* grad3D;
	initializeGradData(grad3D, imageHeight, imageLength, imageWidth, p);
	//==============================================================================
	// Set the ABSContainer
	ABSContainer* dta3D;
	initializeABSContainer(dta3D, imageHeight, imageLength, imageWidth, p);
	//==============================================================================
	// Shorten param names
	dataType foreground = (*dta3D).imgForeground, background = (*dta3D).imgBackground;
	//==============================================================================
	// Estimate the initial segmentation and store to pcaPriorShape - 0th step?
	(*dta3D).est_start = false;
	//==============================================================================
	// Copy the imageData to the (*dta3D).dta - used inside the segmentation fn's
	copyDataPointer(imageData.imageDataPtr, (*dta3D).dta, imageHeight, imageLength, imageWidth);
	//==============================================================================
	// Run the gradient coeffiecient fn
	gradientCoefficients(dta3D, grad3D, imageHeight, imageLength, imageWidth, p);
	//==============================================================================
	// Run the initial segmentation fn
	segmentation3D_fn(dta3D, atls3D, (*dta3D).imgBackground, p, imageHeight, imageLength, imageWidth);
	//==============================================================================
	// Estimate initial segmentation - store as the 0th step priorShape
	// Copy the zeroth or initial segmentation
	copyDataPointer((*dta3D).dta_u, tmpDataStepHolder->d_0th_step, imageHeight + 2 * p + 1, imageLength + 2 * p + 1, imageWidth + 2 * p + 1);
	//==============================================================================
	// Get stored PCA results - mean shape, eigenvectors, eigenvalues, from the pcaParam
	//==============================================================================
	// Read the data for PCA from the files
	readPCAData(pcaParam, pcaParam->pcaFolder, pcaParam->meanDistFile, pcaParam->meanFile, pcaParam->pcaResultFile);
	// Calculate the dist. map for the mean shape
	fastSweepingFunction_3D((*pcaParam).dtaMeanDist, (*pcaParam).dtaMean, imageLength, imageWidth, imageHeight, 1, 50000, foreground);
	ReDistFn((*pcaParam).dtaMeanDist, imageHeight, imageLength, imageWidth, foreground, background);
	//==============================================================================
	// Set PCA paramaters and methods
	dataType total = 0.0, eigsum = 0.0;
	int k_comps = 0;;
	for (int i = 0; i < (*pcaParam).princomp; i++)
	{
		total += fabsf((*pcaParam).eigenvalues[i]);
	}
	for (int i = 0; i < (*pcaParam).princomp; i++)
	{
		eigsum += fabsf((*pcaParam).eigenvalues[i]) / total;
		k_comps++;
		if (eigsum >= pcaParam->threshold_comp)
		{
			break;
		}
	}
	// Set no. of components
	(*pcaParam).princomp = k_comps;
	//==============================================================================
	// Initialize PCA Prior Segmentation shape to distance map for Mean shape
	(*dta3D).priorShape = (dataType**)malloc(mem_alloc_h);
	(*dta3D).estShape = (dataType**)malloc(mem_alloc_h);
	for (size_t i = 0; i < (imageHeight + 2 * p + 1); i++)
	{
		(*dta3D).priorShape[i] = (dataType*)malloc(mem_alloc_lw);
		(*dta3D).estShape[i] = (dataType*)malloc(mem_alloc_lw);
		for (size_t j = 0; j < dimLW; j++)
		{
			(*dta3D).priorShape[i][j] = foreground;
			(*dta3D).estShape[i][j] = foreground;
		}
	}
	//==============================================================================
	// Check if the initial segmentation is similar to the estimated shape
	copyDataPointer((*dta3D).dta_u, (*dta3D).dta_tmp, ((*dta3D).height), ((*dta3D).length), ((*dta3D).width));
	if (estmateSegmentation(dta3D, pcaParam, (*dta3D).dta_tmp, (*dta3D).height, (*dta3D).length, (*dta3D).width, p, foreground, background, tol_s, h3, regParams, optMethod) < tol_e)
	{
		printf("\nEstimated shape found very similar to current segmentation at step: %d\n", 0);
	}
	//==============================================================================
	// Section to read precise result if any, it's dist. fn. This will can be used to compare with final segmentation results
	//==============================================================================
	// Initialize the mass_diff to a large negative value;
	dataType mass_diff = INT_MIN;
	// Loope through the iterations until either reach max no. or tolerance as the stopping condtion
	for (size_t w = 1; w < max_iters; w++)
	{
		//==============================================================================
		// Call the segmentation_ap to calculate the ap eq. coefficients
		segmentation3D_Ap_coef(
			dta3D, grad3D, atls3D,
			imageHeight, imageLength, imageWidth,
			(*dta3D).priorShape,
			lambda, eta, zeta, h3, p,
			epsilon, tau, mu1, mu2, d1, d2
		);
		//==============================================================================
		// Call the gmcf_atlass fn.
		gmcf3D_atlas(
			dta3D, grad3D, atls3D,
			imageHeight, imageLength, imageWidth,
			p
		);
		//==============================================================================
		// Find the error between current segmentation and precise didt. fn's if we have the precise
		//==============================================================================
		/*
		* Call the stop segment function
		*/
		bool stopSegmentation = stop_segment3D(dta3D, tmpDataStepHolder, p, h3, &w, &lambda, &zeta, &eta, &mass_diff, &max_iters, tol_s, tol_e, pcaParam, optMethod, regParams);
		if (stopSegmentation) {
			//==============================================================================
			// Can cacl. the diff between segmentation result and precise result if we have
			//==============================================================================
			// Print the mass_diff
			printf("The final error is %.12e\n", mass_diff);
			//==============================================================================
			// Clean up all pointers created
			// Free dta3d
			// Free atlas3d
			// Free grad3D
			freeAtlasABSContainerGradData(dta3D, atls3D, grad3D, imageHeight, p);
			//==============================================================================
			// Break out of the loop
			break;
			//==============================================================================

		}
	}
}
//==============================================================================
void freeAtlasABSContainerGradData(ABSContainer* dta3D, AtlasData * atls3D, GradData* grad3D, size_t imageHeight, size_t reflexLength)
	{
		size_t k, height = imageHeight + 2 * reflexLength + 1;
		for (k = 0; k < height; k++)
		{
			//==============================================================================
			// dta3D
			free((*dta3D).dta_h[k]);
			free((*dta3D).dta_l[k]);
			free((*dta3D).dta_w[k]);

			free((*dta3D).dta_v[k]);
			free((*dta3D).dta_u[k]);
			free((*dta3D).dta_u0[k]);
			free((*dta3D).dta_uq[k]);
			free((*dta3D).dta_tmp[k]);
			free((*dta3D).dta_tmp2[k]);
			free((*dta3D).dta_tmp3[k]);

			free((*dta3D).priorShape[k]);
			free((*dta3D).estShape[k]);
			//==============================================================================
			// atlas3D
			free((*atls3D).bp_e[k]);
			free((*atls3D).bp_w[k]);
			free((*atls3D).bp_s[k]);
			free((*atls3D).bp_n[k]);
			free((*atls3D).bp_b[k]);
			free((*atls3D).bp_t[k]);
			free((*atls3D).cp[k]);
			free((*atls3D).dp[k]);
			free((*atls3D).gd[k]);
			//==============================================================================
			// grad3D
			free((*grad3D).ae[k]);
			free((*grad3D).aw[k]);
			free((*grad3D).an[k]);
			free((*grad3D).as[k]);
			free((*grad3D).at[k]);
			free((*grad3D).ab[k]);

			free((*grad3D).ge[k]);
			free((*grad3D).gw[k]);
			free((*grad3D).gn[k]);
			free((*grad3D).gs[k]);
			free((*grad3D).gt[k]);
			free((*grad3D).gb[k]);
			free((*grad3D).gd[k]);

			free((*grad3D).ap[k]);
			//==============================================================================
		}
		//==============================================================================
		// Free atlas
		free((*atls3D).bp_e);
		free((*atls3D).bp_w);
		free((*atls3D).bp_s);
		free((*atls3D).bp_n);
		free((*atls3D).bp_b);
		free((*atls3D).bp_t);

		free((*atls3D).cp);
		free((*atls3D).dp);

		free((*atls3D).gd);

		// Set to NULL

		(*atls3D).bp_e = NULL;
		(*atls3D).bp_w = NULL;
		(*atls3D).bp_s = NULL;
		(*atls3D).bp_n = NULL;
		(*atls3D).bp_b = NULL;
		(*atls3D).bp_t = NULL;

		(*atls3D).cp = NULL;
		(*atls3D).dp = NULL;

		(*atls3D).gd = NULL;

		// Free Structs
		free(atls3D->bp_e);
		free(atls3D->bp_w);
		free(atls3D->bp_s);
		free(atls3D->bp_n);
		free(atls3D->bp_b);
		free(atls3D->bp_t);

		free(atls3D->cp);
		free(atls3D->dp);

		free(atls3D->gd);;

		free(atls3D);

		atls3D = NULL;
		//==============================================================================
		// Free grad3D
		free((*grad3D).ae);
		free((*grad3D).aw);
		free((*grad3D).an);
		free((*grad3D).as);
		free((*grad3D).at);
		free((*grad3D).ab);

		free((*grad3D).ge);
		free((*grad3D).gw);
		free((*grad3D).gn);
		free((*grad3D).gs);
		free((*grad3D).gt);
		free((*grad3D).gb);
		free((*grad3D).gd);

		free((*grad3D).ap);

		// Set to Null

		(*grad3D).ae = NULL;
		(*grad3D).aw = NULL;
		(*grad3D).an = NULL;
		(*grad3D).as = NULL;
		(*grad3D).at = NULL;
		(*grad3D).ab = NULL;

		(*grad3D).ge = NULL;
		(*grad3D).gw = NULL;
		(*grad3D).gn = NULL;
		(*grad3D).gs = NULL;
		(*grad3D).gt = NULL;
		(*grad3D).gb = NULL;
		(*grad3D).gd = NULL;

		(*grad3D).ap = NULL;
		
		free(grad3D);
		grad3D = NULL;

		//==============================================================================
		// Free dta3D
		free((*dta3D).dta_h);
		free((*dta3D).dta_l);
		free((*dta3D).dta_w);

		free((*dta3D).dta_v);
		free((*dta3D).dta_u);
		free((*dta3D).dta_u0);
		free((*dta3D).dta_uq);
		free((*dta3D).dta_tmp);
		free((*dta3D).dta_tmp2);
		free((*dta3D).dta_tmp3);

		free((*dta3D).priorShape);
		free((*dta3D).estShape);

		// Set to NULL
		(*dta3D).dta_h = NULL;
		(*dta3D).dta_l = NULL;
		(*dta3D).dta_w = NULL;

		(*dta3D).dta_v = NULL;
		(*dta3D).dta_u = NULL;
		(*dta3D).dta_u0 = NULL;
		(*dta3D).dta_uq = NULL;
		(*dta3D).dta_tmp = NULL;
		(*dta3D).dta_tmp2 = NULL;
		(*dta3D).dta_tmp3 = NULL;

		(*dta3D).priorShape = NULL;
		(*dta3D).estShape = NULL;
		
		free(dta3D);
		dta3D = NULL;
		//==============================================================================
	}
//==============================================================================