#include "front_propagation.h"

typedef struct {
	dataType tau; // used in the descent gradient
	size_t max_iteration; // maximal iteration to stop descent gradient
	dataType tolerance; //minimal distance to stop
} Path_Parameters;

bool shortestPath2d(Image_Data2D actionMapStr, Point2D* seedPoints, std::vector<Point2D>& path_points, Path_Parameters parameters);