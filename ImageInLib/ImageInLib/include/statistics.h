#pragma once
#include "shapeAnalysis_params.h"

void calcMeanShape(dataType **, Shapes *shapes, size_t height, size_t length, size_t width, size_t numShapes, shapeAnalysisParameters params);
void estimateSimilarShape(Shapes * inputShapes, dataType **, shapeAnalysisParameters shapeParam, dataType pcaThreshold, size_t height, size_t length, size_t width, size_t numShapes);
