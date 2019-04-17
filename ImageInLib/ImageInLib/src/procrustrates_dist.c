//==============================================================================
#include <stdio.h>
#include <stdlib.h>
//==============================================================================
#include "procrustrates_dist.h"
//==============================================================================
void calc_mean(dataType **dtaMean, dataType ** dtaInput, size_t height, size_t length, size_t width, size_t numShapes);
void multiplication(dataType **arr1, dataType **arr2, dataType **arr3, const size_t m, const size_t n, const size_t n1);
void transpose(dataType **tposed, dataType **entry, const size_t m, const size_t n);
void copyShapes(dataType * eigshape, dataType **shape, size_t height, size_t length, size_t width);
void shapeEstimate(dataType ** dtaMeanShape, dataType ** shape, dataType * eigenvalues, dataType ** eigenvectors, const size_t princomp, const size_t height, const size_t length, const size_t width);
//==============================================================================
void genProcMeanShape(dataType ** dtaMnShp, Shapes *shapes, size_t height, size_t length, size_t width, size_t numShapes, shape_Analysis_Parameters params)
{
	size_t k, i, j, xd;
	char st[160];
	for (i = 0; i < 160; i++)
	{
		st[i] = '=';
	}
	st[159] = '\0';
	// Initialize the dtaMeanShape and initial estimate shape
	dataType ** initEstimate = (dataType **)malloc(sizeof(dataType*)*height);
	for (k = 0; k < height; k++)
	{
		initEstimate[k] = (dataType*)malloc(sizeof(dataType)*(length*width));
	}
	//==============================================================================
	// Choose and initial estimate of the mean shape
	for (k = 0; k < height; k++)
	{
		for (i = 0; i < length; i++)
		{
			for (j = 0; j < width; j++)
			{
				// 1D Conversion of row and column
				xd = x_new(i, j, length);
				initEstimate[k][xd] = shapes[params.initEstimateShape].shp[k][xd]; // 1st shape as the estimate initial
			}
		}
	}
	//==============================================================================
	// Step two - Register/Align the eigen shapes with the selected mean shape
	int z;
	dataType error, tol_proc = 1.0e-02;
	z = 0;
	do
	{
		//==============================================================================
		z++;
		//==============================================================================
		// Register the shapes with the selected mean shape
		printf("%s\n", st);
		printf("Iteration %d\n", z);
		printf("%s\n\n", st);
		//==============================================================================
		for (k = 0; k < numShapes; k++)
		{
			printf("Shape %zd\n", k + 1);
			run_registration(shapes[k].shp, initEstimate, shapes[k].shp, height, length, width, params.regParams, params.gdescentMethod);
			// Calculate the Mean shape
			calc_mean(dtaMnShp, shapes[k].shp, height, length, width, numShapes);
			printf("%s\n", st);
			printf("\n");
		}
		//==============================================================================
		error = energyFunction(dtaMnShp, initEstimate, height, length, width, params.regParams.h);
		printf("Error is %e\n", error);
		//==============================================================================
		// Update the new estimated mean
		for (k = 0; k < height; k++)
		{
			for (i = 0; i < length; i++)
			{
				for (j = 0; j < width; j++)
				{
					// 1D Conversion of row and column
					xd = x_new(i, j, length);

					initEstimate[k][xd] = dtaMnShp[k][xd];
				}
			}
		}
		// Free memory
		//==============================================================================
	} while ((error > tol_proc) && (z < params.meanCalcSteps));
	printf("The number of iterations to calc. mean shape is %d\n", z);
	printf("Error is %e\n", error);
	// Free memory
	for (k = 0; k < height; k++)
	{
		free(initEstimate[k]);
	}
	free(initEstimate);
	initEstimate = NULL;
}
void pca_analysis(dataType ** dtaMeanShape, dataType *** eigvectors, dataType ** eigvalues, size_t * princomp, Shapes *shapes, size_t height, size_t length, size_t width, size_t numShapes, dataType pca_Threshold)
{
	size_t k, i, j, l, xd, xyd, dim3D = height * length * width;
	// Initialize pointer functions
	dataType **matrix(), *vector();
	// Calculate the center shapes according to mean shape
	// Create saved pointers for the alligned shapes
	dataType ** s_Shapes = (dataType**)malloc(sizeof(dataType*) * numShapes); // S_Transpoe - n by D
	for (i = 0; i < numShapes; i++)
	{
		s_Shapes[i] = malloc(sizeof(dataType) * dim3D);
	}
	//==============================================================================
	//Create feature vector for shapes - D by 1 dimension
	Eigen_Shapes *eigshapes = (Eigen_Shapes*)malloc(sizeof(Eigen_Shapes)*numShapes);
	// Alloc. and copy
	for (i = 0; i < numShapes; i++)
	{
		eigshapes[i].eigenShape = malloc((dim3D) * sizeof(dataType**));
		// Copy values
		copyShapes(eigshapes[i].eigenShape, shapes[i].shp, height, length, width);
	}
	//==============================================================================
	// Fiil in the data from shapes - S transpose
	for (l = 0; l < numShapes; l++)
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

					// center means
					s_Shapes[l][xyd] = eigshapes[l].eigenShape[xyd]/* - dtaMeanShape[k][xd]*/;
				}
			}
		}
	}
	// Free
	free(eigshapes->eigenShape);
	free(eigshapes);
	eigshapes = NULL;
	// Transposed matrix
	dataType ** transp = (dataType**)malloc(sizeof(dataType*) * dim3D); // S - D by n
	for (i = 0; i < dim3D; i++)
	{
		transp[i] = malloc(sizeof(dataType) * numShapes);
	}
	transpose(transp, s_Shapes, numShapes, dim3D);
	// Free
	for (i = 0; i < numShapes; i++)
	{
		free(s_Shapes[i]);
	}
	free(s_Shapes);
	s_Shapes = NULL;
	//==============================================================================
	// Calc the covariance matrix
	dataType **data = matrix(numShapes, dim3D);
	for (i = 1; i <= dim3D; i++)
	{
		for (j = 1; j <= numShapes; j++)
		{
			data[i][j] = transp[i - 1][j - 1];
		}
	}
	// Eigenvectors and Eigenvalues reduction, calculation
	// Parameters
	//==============================================================================
	dataType **covmat_L = matrix(numShapes, numShapes);
	// Covariance calc.
	covcol(data, dim3D, numShapes, covmat_L);
	// Free
	free_matrix(data, dim3D, numShapes);
	//==============================================================================
	dataType *eig_values, *interim;
	eig_values = vector(numShapes);     /* Storage alloc. for vector of eigenvalues */
	interim = vector(numShapes);    /* Storage alloc. for 'intermediate' vector */
	//==============================================================================
	triDecomp(covmat_L, numShapes, eig_values, interim);  /* Triangular decomposition */
	triDian(eig_values, interim, numShapes, covmat_L);   /* Reduction of sym. trid. matrix */
	//==============================================================================
	free_vector(interim, numShapes);
	//==============================================================================
	// Sort the eigen values/vectors in descending order ?
	eigenSort(eig_values, covmat_L, numShapes);
	printf("\nSorted Eigenvalues:\n");
	for (j = numShapes; j >= 1; j--)
	{
		printf("%10.2f\n", eig_values[j]);
	}
	//==============================================================================
	// Put a check to view the corresponding eigenvectors
	printf("\nSorted Eigenvectors:\n");
	for (j = 1; j <= numShapes; j++) {
		for (i = 1; i <= numShapes; i++)
		{
			printf("%12.6f", covmat_L[j][numShapes - i + 1]);
		}
		printf("\n");
	}
	//==============================================================================
	// Selct only significant eigenvalues, sum > 90%
	(*princomp) = 0;
	dataType eigsum = 0.0, total = 0;
	for (i = numShapes; i >= 1; i--)
	{
		total += fabs(eig_values[i]);
	}
	for (i = 1; i <= numShapes; i++)
	{
		eigsum += fabs(eig_values[i]) / total;
		(*princomp)++;
		if (eigsum >= pca_Threshold)
		{
			break;
		}
	}
	//==============================================================================
	printf("\nSelected Eigenvalues:\n");
	for (j = 1; j <= (*princomp); j++)
	{
		printf("%10.2f\n", eig_values[j]);
	}
	//==============================================================================
	// Put a check to view the corresponding eigenvectors
	printf("\nSelected Eigenvectors:\n");
	for (j = 1; j <= numShapes; j++) {
		for (i = numShapes; i > numShapes - (*princomp); i--)
		{
			printf("%14.8f", covmat_L[j][numShapes - i + 1]);
		}
		printf("\n");
	}
	printf("\nWe have selected %zd principal components\n", (*princomp));
	//==============================================================================
	// Calculate the Eigenvalues and eigenvectos of S*S_transpose
	// 1. Align the S - D by n matrix
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
				for (l = 0; l < numShapes; l++)
				{
					// center means
					transp[xyd][l] = transp[xyd][l] - dtaMeanShape[k][xd];
				}
			}
		}
	}
	// 2. Matrix Multiplication to save the eigenvector
	for (i = 0; i < dim3D; i++)
	{
		(*eigvectors)[i] = (dataType*)malloc(sizeof(dataType) * (*princomp));
	}
	(*eigvalues) = (dataType*)malloc(sizeof(dataType) * (*princomp));
	for (j = 0; j < (*princomp); j++)
	{
		(*eigvalues)[j] = eig_values[j + 1];
	}
	free_vector(eig_values, numShapes);
	dataType **selected_eigvectors = (dataType**)malloc(sizeof(dataType*) * numShapes);
	for (i = 0; i < numShapes; i++)
	{
		selected_eigvectors[i] = (dataType*)malloc(sizeof(dataType) * (*princomp));
	}
	for (j = 1; j <= numShapes; j++) {
		for (i = numShapes; i > numShapes - (*princomp); i--)
		{
			selected_eigvectors[j - 1][numShapes - i] = covmat_L[j][numShapes - i + 1];
		}
	}
	free_matrix(covmat_L, numShapes, numShapes); // Free pointer no longer being used
	multiplication(transp, selected_eigvectors, (*eigvectors), dim3D, numShapes, (*princomp));
	// Free
	for (i = 0; i < dim3D; i++)
	{
		free(transp[i]);
	}
	free(transp);
	transp = NULL;
	//==============================================================================
	// Estimate shape
	// Testing
	//estimatedShape(dtaMeanShape, shapes[1].shp, (*eigvalues), (*eigvectors), (*princomp), height, length, width);
	//==============================================================================
	for (i = 0; i < numShapes; i++)
	{
		free(selected_eigvectors[i]);
	}
	free(selected_eigvectors);
	selected_eigvectors = NULL;
	//==============================================================================
}
void estimateShape(Shapes * atlasShapes, dataType ** shapeToEstimate, shape_Analysis_Parameters shapeParam, size_t height, size_t length, size_t width, size_t numShapes, dataType pca_Threshold)
{
	size_t k;
	// Init mean shape pointer
	// Initialize the dtaMeanShape and initial estimate shape
	dataType ** meanShape = (dataType **)malloc(sizeof(dataType*)*height);
	for (k = 0; k < height; k++)
	{
		meanShape[k] = (dataType*)malloc(sizeof(dataType)*(length*width));
	}
	//==============================================================================
	initialize3dArrayD(meanShape, length, width, height, 0);
	//==============================================================================
	// Calculate the mean shape
	genProcMeanShape(meanShape, atlasShapes, height, length, width, numShapes, shapeParam);
	//==============================================================================
	size_t dim3D = height * length*width;
	// PCA Params
	PCA_Params *pcaParam = (PCA_Params*)malloc(sizeof(PCA_Params));
	(*pcaParam).eigenvectors = (dataType**)malloc(sizeof(dataType*) * dim3D);
	(*pcaParam).eigenvalues = NULL;
	// Use Mean shape to calculate the PCA Variables - Eigenvectors, Eigenvalues, Number of Proncipal components
	pca_analysis(meanShape, &(*pcaParam).eigenvectors, &(*pcaParam).eigenvalues, &(*pcaParam).princomp, atlasShapes, height, length, width, numShapes, pca_Threshold);
	//==============================================================================
	// Estimate the Given Shape using the evaluated Pca Parameters
	//==============================================================================
	shapeEstimate(meanShape, shapeToEstimate, (*pcaParam).eigenvalues, (*pcaParam).eigenvectors, (*pcaParam).princomp, height, length, width);
	//==============================================================================
	free((*pcaParam).eigenvectors);
	free(pcaParam);
	//==============================================================================
	for (k = 0; k < height; k++)
	{
		free(meanShape[k]);
	}
	free(meanShape);
	meanShape = NULL;
	//==============================================================================
}
//==============================================================================
void calc_mean(dataType ** dtaMean, dataType ** dtaInput, size_t height, size_t length, size_t width, size_t numShapes)
{
	size_t k, i, j, xd;

	for (k = 0; k < height; k++)
	{
		for (i = 0; i < length; i++)
		{
			for (j = 0; j < width; j++)
			{
				// 2D to 1D representation for i, j
				xd = x_new(i, j, length);

				dtaMean[k][xd] += (dtaInput[k][xd]) / numShapes;
			}
		}
	}
}
//==============================================================================
void transpose(dataType **tposed, dataType **entry, const size_t m, const size_t n)
{
	int j, i;
	// Interchange n with m, j with i
	for (j = 0; j < n; ++j)
	{
		for (i = 0; i < m; ++i)
		{
			tposed[j][i] = entry[i][j];
		}
	}
} // End Matrix Transpose Function
void copyShapes(dataType * eigshape, dataType **shape, size_t height, size_t length, size_t width)
{
	size_t k, i, j, xd, xyd;

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

				eigshape[xyd] = shape[k][xd];
			}
		}
	}
}
void shapeEstimate(dataType ** dtaMeanShape, dataType ** shape, dataType * eigenvalues, dataType ** eigenvectors, const size_t princomp, const size_t height, const size_t length, const size_t width)
{
	size_t k, i, j, xd, xyd;
	//==============================================================================
	size_t dim3D = height * length * width;
	//==============================================================================
	// Create a diagonal eigenvalues matrix inverse
	dataType **diag_eigvalues_inv = (dataType**)malloc(sizeof(dataType*) * princomp); // K by k
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
	dataType ** trans_eigvectors = (dataType**)malloc(sizeof(dataType*) * princomp); // K by D
	for (i = 0; i < princomp; i++)
	{
		trans_eigvectors[i] = (dataType*)malloc(sizeof(dataType) * dim3D);
	}
	transpose(trans_eigvectors, eigenvectors, dim3D, princomp);
	//==============================================================================
	// Difference btn shape and mean shape - D by 1
	dataType * shp_diff = (dataType*)malloc(sizeof(dataType) * dim3D); // D by 1
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
	dataType * shp_projection = (dataType*)malloc(sizeof(dataType) * princomp); // k by 1
	dataType * tmp_results = (dataType*)malloc(sizeof(dataType) * princomp); // k by 1
	// Multiplication
	for (i = 0; i < princomp; i++)
	{
		shp_projection[i] = 0;
		for (j = 0; j < dim3D; j++)
		{
			shp_projection[i] += trans_eigvectors[i][j] * shp_diff[j];
		}
		// Is this Thresholding applicable to general data?
		if (shp_projection[i] > 3)
		{
			shp_projection[i] = 3;
		}
		else if (shp_projection[i] < -3)
		{
			shp_projection[i] = -3;
		}
	}
	//==============================================================================
	// Estimated shape parameter - 1 by 1
	dataType energy_param = 0;
	// Shape projection transpose
	dataType ** shp_projection_trans = (dataType**)malloc(sizeof(dataType *) * 1); // 1 by k
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
		for (j = 0; j < princomp; j++)
		{
			tmp_results[i] = diag_eigvalues_inv[i][j] * shp_projection[j];
		}
	}
	//==============================================================================
	for (i = 0; i < princomp; i++)
	{
		energy_param += shp_projection_trans[0][i] * tmp_results[i];
	}
	//==============================================================================
	printf("\nEnergy of Shape Parameter: %.8lf\n", energy_param);
	//==============================================================================
	// Minimization of Energy
	int iteration = 0, steps = 100;
	dataType tolerance = 1.0e-04, error;
	dataType eps = 1000;
	// Initialize pointers#
	dataType * est_shape = (dataType*)malloc(sizeof(dataType*) * dim3D); // S - D by 1 components
	for (i = 0; i < dim3D; i++)
	{
		est_shape[i] = 0;
	}
	//==============================================================================
	do
	{
		//==============================================================================
		// print shp projection
		printf("\nShape Projection vector: ");
		for (j = 0; j < princomp; j++) {
			printf("%.4f\n", shp_projection[j]);
		}
		//==============================================================================
		printf("Iteration: %d\n", iteration);
		//==============================================================================
		// Calc. Estimated U_k*a_k
		for (i = 0; i < dim3D; i++)
		{
			for (j = 0; j < princomp; j++)
			{
				// Estimated shape
				est_shape[i] += shp_projection[j] * sqrt(eigenvalues[j]) * eigenvectors[i][j];
			}
		}
		//==============================================================================
		// Error
		error = 0;
		for (i = 0; i < princomp; i++)
		{
			error += shp_projection_trans[0][i] * tmp_results[i];
		}
		//error = energyFunction(dtaMeanShape, prior_Shape, height, length, width, h);
		printf("Error is %.8lf\n", error);
		//==============================================================================
		for (i = 0; i < princomp; i++)
		{
			shp_projection[i] = shp_projection[i] - eps * tmp_results[i];
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
			for (j = 0; j < princomp; j++)
			{
				tmp_results[i] = diag_eigvalues_inv[i][j] * shp_projection[j];
			}
		}
		//==============================================================================
		iteration++;
		//==============================================================================
	} while ((error > tolerance) && (iteration < steps));
	printf("\nThe number of iterations to calc. prior shape is %d\n", iteration);
	printf("\nError during calc.: %.8lf\n", error);
	//==============================================================================
	// Rewrite the to old estmate shape
	// Rescale the data?
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
				shape[k][xd] = dtaMeanShape[k][xd] + est_shape[xyd];
			}
		}
	}
	//==============================================================================
	// Free memory
	free(est_shape);
	est_shape = NULL;
	//==============================================================================
	// Free Memory
	for (i = 0; i < princomp; i++)
	{
		free(diag_eigvalues_inv[i]);
		free(trans_eigvectors[i]);
	}
	free(diag_eigvalues_inv);
	free(trans_eigvectors);
	free(tmp_results);
	free(shp_projection);

	tmp_results = NULL;
	shp_projection = NULL;
	diag_eigvalues_inv = NULL;
	trans_eigvectors = NULL;
	for (i = 0; i < 1; i++)
	{
		free(shp_projection_trans[i]);
	}
	free(shp_projection_trans);
	free(shp_diff);
	free(shp_projection);
	//==============================================================================
}
//==============================================================================
// Perform two matrices multiplication Function
void multiplication(dataType **arr1, dataType **arr2, dataType **arr3, const size_t  m, const size_t  n, const size_t n1)
{
	int i, j, k;
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n1; j++)
		{
			arr3[i][j] = 0;
			for (k = 0; k < n; k++)
			{
				arr3[i][j] = arr3[i][j] + arr1[i][k] * arr2[k][j];
			}
		}
	}
} // End Matrices muplication Function
//==============================================================================
dataType procDist(dataType ** dta1, dataType ** dta2, size_t height, size_t length, size_t width)
{
	size_t k, i, j, xd;
	dataType procTmp = 0.0;
	dataType x = 0.0, y = 0.0, z = 0.0;
	for (k = 0; k < height; k++)
	{
		for (i = 0; i < length; i++)
		{
			for (j = 0; j < width; j++)
			{
				// 2D to 1D representation for i, j
				xd = x_new(i, j, length);

				procTmp += (dta1[k][xd] - dta2[k][xd])*(dta1[k][xd] - dta2[k][xd]);
			}
		}
	}
	return sqrt(procTmp);
}
//==============================================================================