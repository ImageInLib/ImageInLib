#include"../src/data_storage.h"
#include "interpolations.h"
#include "imageInterpolation.h"
#include<math.h>

//3D images

Point3D getImageCoordFromRealCoord3D(Point3D image_point, Point3D image_origin, VoxelSpacing image_spacing, OrientationMatrix orientation) {
    Point3D resultPoint;
    resultPoint.x = image_origin.x + image_spacing.sx * orientation.v1.x * image_point.x;
    resultPoint.y = image_origin.y + image_spacing.sy * orientation.v2.y * image_point.y;
    resultPoint.z = image_origin.z + image_spacing.sz * orientation.v3.z * image_point.z;
    return resultPoint;
}

Point3D getRealCoordFomImageCoord3D(Point3D real_point, Point3D real_origin, VoxelSpacing real_spacing, OrientationMatrix orientation) {
    Point3D resultPoint;
    resultPoint.x = (real_point.x - real_origin.x) / real_spacing.sx * orientation.v1.x;
    resultPoint.y = (real_point.y - real_origin.y) / real_spacing.sy * orientation.v2.y;
    resultPoint.z = (real_point.z - real_origin.z) / real_spacing.sz * orientation.v3.z;
    return resultPoint;
}

//=======================================================
//2D images

Point2D getRealCoordFromImageCoord2D(Point2D image_point, Point2D real_origin, PixelSpacing real_spacing, OrientationMatrix2D orientation) {
    Point2D result_point;
    result_point.x = real_origin.x + real_spacing.sx * orientation.v1.x * image_point.x;
    result_point.y = real_origin.y + real_spacing.sy * orientation.v2.y * image_point.y;
    return result_point;
}

Point2D getImageCoordFromRealCoord2D(Point2D real_point, Point2D real_origin, PixelSpacing real_spacing, OrientationMatrix2D orientation) {
    Point2D result_point;
    result_point.x = (real_point.x - real_origin.x) / real_spacing.sx * orientation.v1.x;
    result_point.y = (real_point.y - real_origin.y) / real_spacing.sy * orientation.v2.y;
    return result_point;
}

PointNeighbors2D getPointNeighbors2D(Point2D point, PixelSpacing spacing) {

    PointNeighbors2D neighbors;

    int i_floor = floor(point.x / spacing.sx), i_ceil = ceil(point.x / spacing.sx);
    int j_floor = floor(point.y / spacing.sy), j_ceil = ceil(point.y / spacing.sy);

    neighbors.top_left.x = i_floor * spacing.sx; neighbors.top_left.y = j_floor * spacing.sy;
    neighbors.bottom_left.x = i_floor * spacing.sx; neighbors.bottom_left.y = j_ceil * spacing.sy;
    neighbors.top_right.x = i_ceil * spacing.sx; neighbors.top_right.y = j_floor * spacing.sy;
    neighbors.bottom_right.x = i_ceil * spacing.sx; neighbors.bottom_left.y = j_ceil * spacing.sy;

    return neighbors;
}

dataType getDistance2D(Point2D point1, Point2D point2) {
    return sqrt(pow(point1.x - point2.x, 2) + pow(point1.y - point2.y, 2));
}

Point2D getNearestNeighbor2D(Point2D point, PixelSpacing spacing) {

    PointNeighbors2D neighbors;

    //Point2D result_point;

    neighbors = getPointNeighbors2D(point, spacing);

    dataType distance[4];
    distance[0] = getDistance2D(point, neighbors.top_left);
    distance[1] = getDistance2D(point, neighbors.top_right);
    distance[2] = getDistance2D(point, neighbors.bottom_left);
    distance[3] = getDistance2D(point, neighbors.bottom_right);

    dataType min_dist = 1000000000000;
    int k, min_indice;

    for (k = 0; k < 4; k++) {
        if (distance[k] < min_dist) {
            min_dist = distance[k];
            min_indice = k;
        }
    }

    if (min_indice == 0) {
        return neighbors.top_left;
    }
    if (min_indice == 1) {
        return neighbors.top_right;
    }
    if (min_indice == 2) {
        return neighbors.bottom_left;
    }
    if (min_indice == 3) {
        return neighbors.bottom_right;
    }
}

bool interpolateImageFromImageCStoRealWorldCS(Image_Data2D src_image, Image_Data2D interp_image, dataType scale_factor) {

    //The structure interp_image represents the image in real world coordinates system
    //the spacing is the real spacing
    //the origin is the real origin
    //the orientation is the real orientation

    if (src_image.imageDataPtr == NULL || interp_image.imageDataPtr == NULL)
        return false;

    const size_t src_height = src_image.height, src_width = src_image.width;
    const size_t interp_height = interp_image.height, interp_width = interp_image.width;

    int i, j;
    Point2D src_image_point, interp_image_point, src_rw_point, interp_rw_point;

    for (i = 0; i < interp_height; i++) {
        for (j = 0; j < interp_width; j++) {

            //Interpolated image point
            interp_image_point.x = i;
            interp_image_point.y = j;

            //Get the source point correponding to interpolated point
            src_image_point.x = interp_image_point.x * (1. / scale_factor);
            src_image_point.y = interp_image_point.y * (1. / scale_factor);

            //Get src point in real world CS
            src_rw_point = getRealCoordFromImageCoord2D(src_image_point, interp_image.origin, interp_image.spacing, interp_image.orientation);

            //Get the nearest neighbor in real world coordinates system
            //The nearest neighbor is found with real world cordinate system
            Point2D nearest_neighbor = getNearestNeighbor2D(src_rw_point, interp_image.spacing);

            //Get the corresponding point to the nearest neighbor in image coordinates system
            interp_image_point = getImageCoordFromRealCoord2D(nearest_neighbor, interp_image.origin, interp_image.spacing, interp_image.orientation);

            //Set the interpolated value
            int i_int = (int)interp_image_point.x, j_int = (int)interp_image_point.y;
            if (i_int > src_height - 1) {
                i_int = src_height - 1;
            }
            if (j_int > src_width - 1) {
                j_int = src_width - 1;
            }

            interp_image.imageDataPtr[x_new(i, j, interp_height)] = src_image.imageDataPtr[x_new(i_int, j_int, src_height)];
        }
    }
    return true;
}

