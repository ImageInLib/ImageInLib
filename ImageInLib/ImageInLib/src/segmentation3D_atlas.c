#include "segmentation3D_atlas.h"
#include <shape_registration.h>
#include <filter_params.h>
#include <climits>
//==============================================================================
// Function to calclate the gamma parameter - eq 91
dataType chooseGamma(dataType value, dataType dist1, dataType dist2);
//==============================================================================
// Intitialize the atlas data containers
void initializeAtlasData(AtlasData* atls3D, size_t height, size_t length, size_t width, size_t reflexLength);
// Initialize Gradient data containers
void initializeGradData(GradData* grad3D, size_t imageHeight, size_t imageLength, size_t imageWidth, size_t reflexLength);
// Initialize the ABSContainer
void initializeABSContainer(ABSContainer* dta3D, size_t imageHeight, size_t imageLength, size_t imageWidth, size_t reflexLength);
// Estimate shape function by minimization of the energy
void estimatedShapeByEnergy(ABSContainer* dta3D, dataType** dtaMeanShape, dataType** shape, dataType** currentSeg, dataType* eigenvalues, dataType** eigenvectors, int princomp, size_t height, size_t length, size_t width, int step, Affine_Parameter finalResults, dataType bound, dataType eps, int est_iterations, bool parallelize, int NUMTHREADS, dataType h, dataType tolerance);
// Estimate shape function by probability function
void estimateShapeByProbability(ABSContainer* dta3D, dataType** dtaMeanShape, dataType** shape, dataType** currentSeg, dataType* eigenvalues, dataType** eigenvectors, int princomp, size_t height, size_t length, size_t width, int step, Affine_Parameter finalResults, dataType bound, dataType eps, int est_iterations, bool parallelize, int NUMTHREADS, dataType h, dataType tolerance);
//==============================================================================


/* Function implementations*/
//==============================================================================
dataType chooseGamma(dataType value, dataType dist1, dataType dist2)
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
void estimatedShapeByEnergy(ABSContainer* dta3D,dataType** dtaMeanShape, dataType** shape, dataType** currentSeg, dataType* eigenvalues, dataType** eigenvectors, int princomp, size_t height, size_t length, size_t width, int step, Affine_Parameter finalResults, dataType bound, dataType eps, int est_iterations, bool parallelize, int NUMTHREADS, dataType h, dataType tolerance)
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
void estimateShapeByProbability(ABSContainer* dta3D, dataType** dtaMeanShape, dataType** shape, dataType** currentSeg, dataType* eigenvalues, dataType** eigenvectors, int princomp, size_t height, size_t length, size_t width, int step, Affine_Parameter finalResults, dataType bound, dataType eps, int est_iterations, bool parallelize, int NUMTHREADS, dataType h, dataType tolerance)
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
dataType estmateSegmentation(ABSContainer* dta3D, PCAData* pcaParam, dataType** currentSegmentation, size_t height, size_t length, size_t width, size_t reflexLength, dataType foreground, dataType background, dataType tol_s, int step, dataType h)
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
		// ToDo - update the run_registration function to return the optimal registration parameters.
		Affine_Parameter finalResults = run_registration((*pcaParam).dtaMean, (*dta3D).dta_tmp3, regCurrentSegmentation, height, length, width, tol_est, step_size, false);
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
			estimatedShapeByEnergy(dta3D,(*pcaParam).dtaMeanDist, regCurrSegDistMap, currentSegmentation, (*pcaParam).eigenvalues, (*pcaParam).eigenvectors, (*pcaParam).princomp, height, length, width, step, finalResults, bound, eps, est_iterations, false, 2, h, tolerance);
		}
		else if ((*pcaParam).estMethod == MIN_Probability) // Minimizing probability
		{
			bound = 0.25, eps = 1.0;
			estimateShapeByProbability(dta3D ,(*pcaParam).dtaMeanDist, regCurrSegDistMap, currentSegmentation, (*pcaParam).eigenvalues, (*pcaParam).eigenvectors, (*pcaParam).princomp, height, length, width, step, finalResults, bound, eps, est_iterations, false, 2, h, tolerance);
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
void segmentation3D_Ap_coef(ABSContainer* dta3D, GradData* grad3D, AtlasData* atls3D, size_t height, size_t length, size_t width, dataType** priorShape, dataType lambda, dataType eta, dataType zeta, dataType h, size_t p, dataType epsilon, dataType tau, dataType mu1, dataType mu2)
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
				gama = GammaFun(tmp6); // |u - u|

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
				gama = GammaFun(tmp6); // |u - u|

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
				gama = GammaFun(tmp6); // |u - u|

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
				gama = GammaFun(tmp6); // |u - u|

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
				gama = GammaFun(tmp6); // |u - u|

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
				gama = GammaFun(tmp6); // |u - u|

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
// Atlas Segmentation Model Interface
void atlasSegmentationModel(Segmentation_Paramereters segParameters, size_t imageHeight, size_t imageLength, size_t imageWidth)
{
}
