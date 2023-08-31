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

dataType getInterpolatedValueNearestNeighbor3D(Image_Data src_image, Point3D point) {

    int i_floor, i_ceil, j_floor, j_ceil, k_floor, k_ceil;

    if (floor(point.x) < 0) {
        i_floor = 0;
    }
    else {
        if (floor(point.x) > src_image.length - 1) {
            i_floor = src_image.length - 1;
        }
        else {
            i_floor = floor(point.x);
        }
    }

    if (floor(point.y) < 0) {
        j_floor = 0;
    }
    else {
        if (floor(point.y) > src_image.width - 1) {
            j_floor = src_image.width - 1;
        }
        else {
            j_floor = floor(point.y);
        }
    }

    if (floor(point.z) < 0) {
        k_floor = 0;
    }
    else {
        if (floor(point.z) > src_image.height - 1) {
            k_floor = src_image.height - 1;
        }
        else {
            k_floor = floor(point.z);
        }
    }

    if (ceil(point.x) < 0) {
        i_ceil = 0;
    }
    else {
        if (ceil(point.x) > src_image.length - 1) {
            i_ceil = src_image.length - 1;
        }
        else {
            i_ceil = ceil(point.x);
        }
    }

    if (ceil(point.y) < 0) {
        j_ceil = 0;
    }
    else {
        if (ceil(point.y) > src_image.width - 1) {
            j_ceil = src_image.width - 1;
        }
        else {
            j_ceil = ceil(point.y);
        }
    }

    if (ceil(point.z) < 0) {
        k_ceil = 0;
    }
    else {
        if (ceil(point.z) > src_image.height - 1) {
            k_ceil = src_image.height - 1;
        }
        else {
            k_ceil = ceil(point.z);
        }
    }

    dataType x1 = i_floor * src_image.spacing.sx;
    dataType x2 = i_ceil * src_image.spacing.sx;
    dataType y1 = j_floor * src_image.spacing.sy;
    dataType y2 = j_ceil * src_image.spacing.sy;
    dataType z1 = k_floor * src_image.spacing.sz;
    dataType z2 = k_ceil * src_image.spacing.sz;

    //Current point neighbors
    Point3D P000 = { x1, y1, z1 }, P100 = { x2, y1, z1 }, P101 = { x2, y1, z2 };
    Point3D P010 = { x1, y2, z1 }, P011 = { x1, y2, z2 }, P001 = { x1, y1, z2 };
    Point3D P111 = { x2, y2, z2 }, P110 = { x2, y2, z1 };

    double min_dist = 1000000000000;

    int i_int = 0, j_int = 0, k_int = 0;

    double distP000 = getPoint3DDistance(point, P000);
    if (distP000 < min_dist) {
        min_dist = distP000;
        i_int = i_floor;
        j_int = j_floor;
        k_int = k_floor;
    }

    double distP010 = getPoint3DDistance(point, P010);
    if (distP010 < min_dist) {
        min_dist = distP010;
        i_int = i_floor;
        j_int = j_ceil;
        k_int = k_floor;
    }

    double distP100 = getPoint3DDistance(point, P100);
    if (distP100 < min_dist) {
        min_dist = distP100;
        i_int = i_ceil;
        j_int = j_floor;
        k_int = k_floor;
    }

    double distP110 = getPoint3DDistance(point, P110);
    if (distP110 < min_dist) {
        min_dist = distP110;
        i_int = i_ceil;
        j_int = j_ceil;
        k_int = k_floor;
    }

    double distP001 = getPoint3DDistance(point, P001);
    if (distP001 < min_dist) {
        min_dist = distP001;
        i_int = i_floor;
        j_int = j_floor;
        k_int = k_ceil;
    }

    double distP101 = getPoint3DDistance(point, P101);
    if (distP101 < min_dist) {
        min_dist = distP101;
        i_int = i_ceil;
        j_int = j_floor;
        k_int = k_ceil;
    }

    double distP011 = getPoint3DDistance(point, P011);
    if (distP011 < min_dist) {
        min_dist = distP011;
        i_int = i_floor;
        j_int = j_ceil;
        k_int = k_ceil;
    }

    double distP111 = getPoint3DDistance(point, P111);
    if (distP111 < min_dist) {
        min_dist = distP111;
        i_int = i_ceil;
        j_int = j_ceil;
        k_int = k_ceil;
    }

    //set the interpolated value
    return src_image.imageDataPtr[k_int][x_new(i_int, j_int, src_image.length)];

}

dataType getInterpolatedValueTrilinear3D(Image_Data src_image, Point3D point) {

    int i_floor, i_ceil, j_floor, j_ceil, k_floor, k_ceil;

    if (floor(point.x) < 0) {
        i_floor = 0;
    }
    else {
        if (floor(point.x) > src_image.length - 1) {
            i_floor = src_image.length - 1;
        }
        else {
            i_floor = floor(point.x);
        }
    }

    if (floor(point.y) < 0) {
        j_floor = 0;
    }
    else {
        if (floor(point.y) > src_image.width - 1) {
            j_floor = src_image.width - 1;
        }
        else {
            j_floor = floor(point.y);
        }
    }

    if (floor(point.z) < 0) {
        k_floor = 0;
    }
    else {
        if (floor(point.z) > src_image.height - 1) {
            k_floor = src_image.height - 1;
        }
        else {
            k_floor = floor(point.z);
        }
    }

    if (ceil(point.x) < 0) {
        i_ceil = 0;
    }
    else {
        if (ceil(point.x) > src_image.length - 1) {
            i_ceil = src_image.length - 1;
        }
        else {
            i_ceil = ceil(point.x);
        }
    }

    if (ceil(point.y) < 0) {
        j_ceil = 0;
    }
    else {
        if (ceil(point.y) > src_image.width - 1) {
            j_ceil = src_image.width - 1;
        }
        else {
            j_ceil = ceil(point.y);
        }
    }

    if (ceil(point.z) < 0) {
        k_ceil = 0;
    }
    else {
        if (ceil(point.z) > src_image.height - 1) {
            k_ceil = src_image.height - 1;
        }
        else {
            k_ceil = ceil(point.z);
        }
    }

    dataType x1 = (dataType)i_floor;
    dataType x2 = (dataType)i_ceil;
    dataType y1 = (dataType)j_floor;
    dataType y2 = (dataType)j_ceil;
    dataType z1 = (dataType)k_floor;
    dataType z2 = (dataType)k_ceil;

    //Current point neighbors
    Point3D P000 = { x1, y1, z1 }, P100 = { x2, y1, z1 }, P101 = { x2, y1, z2 };
    Point3D P010 = { x1, y2, z1 }, P011 = { x1, y2, z2 }, P001 = { x1, y1, z2 };
    Point3D P111 = { x2, y2, z2 }, P110 = { x2, y2, z1 };

    dataType c000 = src_image.imageDataPtr[k_floor][x_new(i_floor, j_floor, src_image.length)];
    dataType c100 = src_image.imageDataPtr[k_floor][x_new(i_ceil, j_floor, src_image.length)];
    dataType c001 = src_image.imageDataPtr[k_ceil][x_new(i_floor, j_floor, src_image.length)];
    dataType c010 = src_image.imageDataPtr[k_floor][x_new(i_floor, j_ceil, src_image.length)];
    dataType c101 = src_image.imageDataPtr[k_ceil][x_new(i_ceil, j_floor, src_image.length)];
    dataType c110 = src_image.imageDataPtr[k_floor][x_new(i_ceil, j_ceil, src_image.length)];
    dataType c111 = src_image.imageDataPtr[k_ceil][x_new(i_ceil, j_ceil, src_image.length)];
    dataType c011 = src_image.imageDataPtr[k_ceil][x_new(i_floor, j_ceil, src_image.length)];

    //case no interpolation
    if (x1 == x2 && y1 == y2 && z1 == z2) {
        return c000;
    }

    //case linear
    if (x1 == x2 && y1 == y2 && z1 != z2) {
        return linearInterpolation(point.z, z1, z2, c000, c001);
    }
    if (x1 == x2 && y1 != y2 && z1 == z2) {
        return linearInterpolation(point.y, y1, y2, c000, c010);
    }
    if (x1 != x2 && y1 == y2 && z1 == z2) {
        return linearInterpolation(point.x, x1, x2, c000, c100);
    }

    //Case bilinear
    if (x1 == x2 && y1 != y2 && z1 != z2) {
        return bilinearInterpolation(point.y, y1, y2, c000, c011, c101, c110, point.z, z1, z2);
    }
    if (x1 != x2 && y1 != y2 && z1 == z2) {
        return bilinearInterpolation(point.x, x1, x2, c000, c110, c001, c111, point.y, y1, y2);
    }
    if (x1 != x2 && y1 == y2 && z1 != z2) {
        return bilinearInterpolation(point.x, x1, x2, c100, c001, c010, c111, point.z, x1, z2);
    }

    //Case trilinear
    if (x1 != x2 && y1 != y2 && z1 != z2) {
        return trilinearInterpolation(point.x, x1, x2, point.y, y1, y2, c000, c001, c010, c011, c100, c101, c110, c111, point.z, z1, z2);
    }

}

bool imageInterpolation3D(Image_Data src_image, Image_Data dest_image, interpolationMethod method) {

    if (src_image.imageDataPtr == NULL || dest_image.imageDataPtr == NULL)
        return false;

    int i, j, k;

    const size_t src_length = src_image.length, src_width = src_image.width, src_height = src_image.height;
    const size_t dest_length = dest_image.length, dest_width = dest_image.width, dest_height = dest_image.height;

    Point3D current_point;

    for (k = 0; k < dest_height; k++) {
        for (i = 0; i < dest_length; i++) {
            for (j = 0; j < dest_width; j++) {

                current_point.x = (dataType)i;
                current_point.y = (dataType)j;
                current_point.z = (dataType)k;

                //get the corresponding point to the given point in real world coordinates system
                current_point = getRealCoordFromImageCoord3D(current_point, dest_image.origin, dest_image.spacing, dest_image.orientation);

                //get the corresponding point to given point image cs of the source image
                current_point = getImageCoordFromRealCoord3D(current_point, src_image.origin, src_image.spacing, src_image.orientation);

                switch (method)
                {
                case NEAREST_NEIGHBOR:
                    dest_image.imageDataPtr[k][x_new(i, j, dest_length)] = getInterpolatedValueNearestNeighbor3D(src_image, current_point);
                    break;
                case TRILINEAR:
                    dest_image.imageDataPtr[k][x_new(i, j, dest_length)] = getInterpolatedValueTrilinear3D(src_image, current_point);
                    break;
                default:
                    break;
                }
            }
        }
    }
    return true;
}

Statistics getStatisticsPET(Image_Data imageData, Point3D point, dataType radius) {

    size_t i, j, k, nb_point = 0;
    Statistics result;
    Point3D point_pet = getImageCoordFromRealCoord3D(point, imageData.origin, imageData.spacing, imageData.orientation);
    dataType sum = 0.0;
    result.max_data = 0.0; result.min_data = 1000000000000000;

    for (k = 0; k < imageData.height; k++) {
        for (i = 0; i < imageData.length; i++) {
            for (j = 0; j < imageData.width; j++) {
                Point3D current_point = { (dataType)i, (dataType)j, (dataType)k };
                dataType dist = getPoint3DDistance(point_pet, current_point);
                //We are considering just voxel in a ball defined by the radius
                if (dist <= radius) {
                    dataType voxel_value = imageData.imageDataPtr[k][x_new(i, j, imageData.length)];
                    if (voxel_value < result.min_data) {
                        result.min_data = voxel_value;
                    }
                    if (voxel_value > result.max_data) {
                        result.max_data = voxel_value;
                    }
                    sum = sum + voxel_value;
                    nb_point = nb_point + 1;
                }
            }
        }
    }

    result.mean_data = sum / (dataType)nb_point;

    dataType sum_diff = 0.0;
    for (k = 0; k < imageData.height; k++) {
        for (i = 0; i < imageData.length; i++) {
            for (j = 0; j < imageData.width; j++) {
                Point3D current_point = { (dataType)i, (dataType)j, (dataType)k };
                dataType dist = getPoint3DDistance(point_pet, current_point);
                if (dist <= radius) {
                    sum_diff = sum_diff + pow(imageData.imageDataPtr[k][x_new(i, j, imageData.length)] - result.mean_data, 2);
                }
            }
        }
    }

    result.sd_data = sqrt(sum_diff / (dataType)nb_point);

    return result;
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

    if (floor(point.x) < 0) {
        i_floor = 0;
    }
    else {
        if (floor(point.x) > src_image.height - 1) {
            i_floor = src_image.height - 1;
        }
        else {
            i_floor = floor(point.x);
        }
    }

    if (floor(point.y) < 0) {
        j_floor = 0;
    }
    else {
        if (floor(point.y) > src_image.width - 1) {
            j_floor = src_image.width - 1;
        }
        else {
            j_floor = floor(point.y);
        }
    }

    if (ceil(point.x) < 0) {
        i_ceil = 0;
    }
    else {
        if (ceil(point.x) > src_image.height - 1) {
            i_ceil = src_image.height - 1;
        }
        else {
            i_ceil = ceil(point.x);
        }
    }

    if (ceil(point.y) < 0) {
        j_ceil = 0;
    }
    else {
        if (ceil(point.y) > src_image.width - 1) {
            j_ceil = src_image.width - 1;
        }
        else {
            j_ceil = ceil(point.y);
        }
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

    if (floor(point.x) < 0) {
        i_floor = 0;
    }
    else {
        if (floor(point.x) > src_image.height - 1) {
            i_floor = src_image.height - 1;
        }
        else {
            i_floor = floor(point.x);
        }
    }

    if (floor(point.y) < 0) {
        j_floor = 0;
    }
    else {
        if (floor(point.y) > src_image.width - 1) {
            j_floor = src_image.width - 1;
        }
        else {
            j_floor = floor(point.y);
        }
    }

    if (ceil(point.x) < 0) {
        i_ceil = 0;
    }
    else {
        if (ceil(point.x) > src_image.height - 1) {
            i_ceil = src_image.height - 1;
        }
        else {
            i_ceil = ceil(point.x);
        }
    }

    if (ceil(point.y) < 0) {
        j_ceil = 0;
    }
    else {
        if (ceil(point.y) > src_image.width - 1) {
            j_ceil = src_image.width - 1;
        }
        else {
            j_ceil = ceil(point.y);
        }
    }

    dataType x1 = (dataType)i_floor;
    dataType x2 = (dataType)i_ceil;
    dataType y1 = (dataType)j_floor;
    dataType y2 = (dataType)j_ceil;


    //Current point neighbors
    Point2D P11 = { x1, y1 }, P12 = { x2, y1 }, P21 = { x1, y2 }, P22 = { x2, y2 };

    dataType q11 = src_image.imageDataPtr[x_new(i_floor, j_floor, src_image.height)];
    dataType q12 = src_image.imageDataPtr[x_new(i_ceil, j_floor, src_image.height)];
    dataType q21 = src_image.imageDataPtr[x_new(i_floor, j_ceil, src_image.height)];
    dataType q22 = src_image.imageDataPtr[x_new(i_ceil, j_ceil, src_image.height)];

    if (x1 == x2 && y1 == y2) {
        return q11;
    }
    if (x1 == x2 && y1 != y2) {
        return linearInterpolation(point.x, y1, y2, q11, q12);
    }
    if (x1 != x2 && y1 == y2) {
        return linearInterpolation(point.y, x1, x2, q12, q22);
    }
    if (x1 != x2 && y1 != y2) {
        return bilinearInterpolation(point.x, x1, x2, q11, q12, q21, q22, point.y, y1, y2);
    }
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
            case BILINEAR:
                dest_image.imageDataPtr[x_new(i, j, dest_height)] = getInterpolatedValueBilinear2D(src_image, current_point);
                break;
            default:
                break;
            }
        }
    }

    return true;
}
