#include "statistics.h"
#include "procrustrates_dist.h"

void calcMeanShape(dataType **meanShapeOutput, Shapes * inputShapes, size_t height, size_t length, size_t width, size_t numShapes, shape_Analysis_Parameters params)
{
	genProcMeanShape(meanShapeOutput, inputShapes, height, length, width, numShapes, params);
}

void estimateSimilarShape(Shapes * inputShapes, dataType **shapeToEstimate, dataType ** estimatedShape, shape_Analysis_Parameters shapeParam, estimate_Params * estimateParams, dataType pcaThreshold, size_t height, size_t length, size_t width, size_t numShapes)
{
	estimateShape(inputShapes, shapeToEstimate, estimatedShape, shapeParam, estimateParams, height, length, width, numShapes, pcaThreshold);
}