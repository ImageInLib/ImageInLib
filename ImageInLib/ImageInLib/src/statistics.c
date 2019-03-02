#include "statistics.h"
#include "procrustrates_dist.h"

void calcMeanShape(dataType **meanShapeOutput, Shapes * inputShapes, size_t height, size_t length, size_t width, size_t numShapes, shapeAnalysisParameters params)
{
	genProcMeanShape(meanShapeOutput, inputShapes, height, length, width, numShapes, params);
}

void estimateSimilarShape(Shapes * inputShapes, dataType **shapeToEstimate, shapeAnalysisParameters shapeParam, dataType pcaThreshold, size_t height, size_t length, size_t width, size_t numShapes)
{
	estimateShape(inputShapes, shapeToEstimate, shapeParam, height, length, width, numShapes, pcaThreshold);
}