#include "path_extraction.h"

bool shortestPath2d(Image_Data2D actionMapStr, Point2D* seedPoints, std::vector<Point2D>& path_points, Path_Parameters parameters) {

	if (actionMapStr.imageDataPtr == NULL || seedPoints == NULL)
		return false;

	dataType dist_min = 0.0;

	const FiniteVolumeSize2D spacing = { actionMapStr.spacing.sx, actionMapStr.spacing.sy };

	size_t i = (size_t)seedPoints[1].x;
	size_t j = (size_t)seedPoints[1].y;
	dataType x = seedPoints[1].x;
	dataType y = seedPoints[1].y;

	Point2D grad;
	dataType tau = parameters.tau;
	size_t count_iter = 0;

	do {

		getGradient2D(actionMapStr.imageDataPtr, actionMapStr.height, actionMapStr.width, i, j, spacing, &grad);
		dataType gradNorm = sqrt(grad.x * grad.x + grad.y * grad.y);
		if(gradNorm < 1e-6) {
			//Gradient is too small, stop the iteration to avoid division by zero
			break;
		}
		x -= tau * grad.x / gradNorm;
		y -= tau * grad.y / gradNorm;

		Point2D current = { x, y };
		dist_min = getPoint2DDistance(current, seedPoints[0]);
		path_points.push_back(current);

		i = (size_t)round(x);
		j = (size_t)round(y);

		count_iter++;
	} while (dist_min > parameters.tolerance && count_iter < parameters.max_iteration);

	return true;
}
