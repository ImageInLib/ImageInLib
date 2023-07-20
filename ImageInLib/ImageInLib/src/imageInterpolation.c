#include"../src/data_storage.h"
#include "interpolations.h"
#include "imageInterpolation.h"
#include<math.h>

//3D images
Point3D getRealCoordFromImageCoord3D(Point3D image_point, Point3D real_origin, VoxelSpacing real_spacing, OrientationMatrix orientation) {
    Point3D resultPoint;
    resultPoint.x = real_origin.x + real_spacing.sx * orientation.v1.x * image_point.x;
    resultPoint.y = real_origin.y + real_spacing.sy * orientation.v2.y * image_point.y;
    resultPoint.z = real_origin.z + real_spacing.sz * orientation.v3.z * image_point.z;
    return resultPoint;
}

Point3D getImageCoordFromRealCoord3D(Point3D real_point, Point3D real_origin, VoxelSpacing real_spacing, OrientationMatrix orientation) {
    Point3D resultPoint;
    resultPoint.x = (real_point.x - real_origin.x) / (real_spacing.sx * orientation.v1.x);
    resultPoint.y = (real_point.y - real_origin.y) / (real_spacing.sy * orientation.v2.y);
    resultPoint.z = (real_point.z - real_origin.z) / (real_spacing.sz * orientation.v3.z);
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
    result_point.x = (real_point.x - real_origin.x) / (real_spacing.sx * orientation.v1.x);
    result_point.y = (real_point.y - real_origin.y) / (real_spacing.sy * orientation.v2.y);
    return result_point;
}

dataType getInterpolatedValueNearestNeighbor2D(Image_Data2D src_image, Point2D point) {

    int i_floor, i_ceil, j_floor, j_ceil;

    if (floor(point.x / src_image.spacing.sx) > src_image.height - 1) {
        i_floor = src_image.height - 1;
    }
    else {
        i_floor = floor(point.x / src_image.spacing.sx);
    }

    if (floor(point.y / src_image.spacing.sy) > src_image.width - 1) {
        j_floor = src_image.width - 1;
    }
    else {
        j_floor = floor(point.y / src_image.spacing.sy);
    }

    if (ceil(point.x / src_image.spacing.sx) > src_image.height - 1) {
        i_ceil = src_image.height - 1;
    }
    else {
        i_ceil = ceil(point.x / src_image.spacing.sx);
    }

    if (ceil(point.y / src_image.spacing.sy) > src_image.width - 1) {
        j_ceil = src_image.width - 1;
    }
    else {
        j_ceil = ceil(point.y / src_image.spacing.sy);
    }

    dataType x1 = i_floor * src_image.spacing.sx;
    dataType x2 = i_ceil * src_image.spacing.sx;
    dataType y1 = j_floor * src_image.spacing.sy;
    dataType y2 = j_ceil * src_image.spacing.sy;

    //Current point neighbors
    Point2D P11 = { x1, y1 }, P12 = { x2, y1 }, P21 = { x1, y2 }, P22 = { x2, y2 };

    double min_dist = 1000000000000;

    int i_int = 0, j_int = 0;

    double distP11 = getPoint2DDistance(point, P11);
    if (distP11 < min_dist) {
        min_dist = distP11;
        i_int = i_floor;
        j_int = j_floor;
    }

    double distP12 = getPoint2DDistance(point, P12);
    if (distP12 < min_dist) {
        min_dist = distP12;
        i_int = i_floor;
        j_int = j_ceil;
    }

    double distP21 = getPoint2DDistance(point, P21);
    if (distP21 < min_dist) {
        min_dist = distP21;
        i_int = i_ceil;
        j_int = j_floor;
    }

    double distP22 = getPoint2DDistance(point, P22);
    if (distP22 < min_dist) {
        min_dist = distP22;
        i_int = i_ceil;
        j_int = j_ceil;
    }

    //set the interpolated value
    return src_image.imageDataPtr[x_new(i_int, j_int, src_image.height)];
}

dataType getInterpolatedValueBilinear2D(Image_Data2D src_image, Point2D point) {

    int i_floor, i_ceil, j_floor, j_ceil;

    if (floor(point.x / src_image.spacing.sx) > src_image.height - 1) {
        i_floor = src_image.height - 1;
    }
    else {
        i_floor = floor(point.x / src_image.spacing.sx);
    }

    if (floor(point.y / src_image.spacing.sy) > src_image.width - 1) {
        j_floor = src_image.width - 1;
    }
    else {
        j_floor = floor(point.y / src_image.spacing.sy);
    }

    if (ceil(point.x / src_image.spacing.sx) > src_image.height - 1) {
        i_ceil = src_image.height - 1;
    }
    else {
        i_ceil = ceil(point.x / src_image.spacing.sx);
    }

    if (ceil(point.y / src_image.spacing.sy) > src_image.width - 1) {
        j_ceil = src_image.width - 1;
    }
    else {
        j_ceil = ceil(point.y / src_image.spacing.sy);
    }

    dataType x1 = i_floor * src_image.spacing.sx;
    dataType x2 = i_ceil * src_image.spacing.sx;
    dataType y1 = j_floor * src_image.spacing.sy;
    dataType y2 = j_ceil * src_image.spacing.sy;

    //Current point neighbors
    Point2D P11 = { x1, y1 }, P12 = { x2, y1 }, P21 = { x1, y2 }, P22 = { x2, y2 };

    dataType q11 = src_image.imageDataPtr[x_new(i_floor, j_floor, src_image.height)];
    dataType q12 = src_image.imageDataPtr[x_new(i_ceil, j_floor, src_image.height)];
    dataType q21 = src_image.imageDataPtr[x_new(i_floor, j_ceil, src_image.height)];
    dataType q22 = src_image.imageDataPtr[x_new(i_ceil, j_ceil, src_image.height)];

    return bilinearInterpolation(point.x, x1, x2, q11, q12, q21, q22, point.y, y1, y2);

}

bool imageInterpolation2D(Image_Data2D src_image, Image_Data2D dest_image, interpolationMethod method) {
    
    if (src_image.imageDataPtr == NULL || dest_image.imageDataPtr == NULL)
        return false;

    int i, j;

    const size_t src_height = src_image.height, src_width = src_image.width;
    const size_t dest_height = dest_image.height, dest_width = dest_image.width;

    Point2D current_point;

    for (i = 0; i < dest_height; i++) {
        for (j = 0; j < dest_width; j++) {

            current_point.x = (dataType)i;
            current_point.y = (dataType)j;

            //get the corresponding point to the given point in real world coordinates system
            current_point = getRealCoordFromImageCoord2D(current_point, dest_image.origin, dest_image.spacing, dest_image.orientation);

            //get the corresponding point to given point image cs of the source image
            current_point = getImageCoordFromRealCoord2D(current_point, src_image.origin, src_image.spacing, src_image.orientation);

            switch (method)
            {
            case NEAREST_NEIGHBOR:
                dest_image.imageDataPtr[x_new(i, j, dest_height)] = getInterpolatedValueNearestNeighbor2D(src_image, current_point);
                break;
            case BILINEAR :
                dest_image.imageDataPtr[x_new(i, j, dest_height)] = getInterpolatedValueBilinear2D(src_image, current_point);
                break;
            default:
                break;
            }
        }
    }

    return true;
}
