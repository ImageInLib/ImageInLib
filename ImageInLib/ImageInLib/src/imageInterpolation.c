#include"../src/data_storage.h"
#include "interpolations.h"
#include "imageInterpolation.h"
#include<math.h>

bool nearestNeighborInterpolation(dataType** originalImage, dataType** newImage, size_t imageLength, size_t imageWidth, size_t imageHeight, dataType originalSpacing, dataType newSpacing)
{
    if (originalImage == NULL || newImage == NULL)
        return false;

    size_t i, j, k, kn, x;
    dataType k_int, k1, k2;

    k_int = 0; kn = 0;

    for (k = 0; k < imageHeight - 1; k++) {
        k1 = originalSpacing * k;
        k2 = originalSpacing * (k + 1);
        do {
            for (i = 0; i < imageLength; i++) {
                for (j = 0; j < imageWidth; j++) {
                    x = x_new(i, j, imageLength);
                    if (k_int == 0) {
                        newImage[kn][x] = originalImage[k][x];
                    }
                    else {
                        if ((k_int - k1) < (k2 - k_int)) {
                            newImage[kn][x] = originalImage[k][x];
                        }
                        else {
                            newImage[kn][x] = originalImage[k + 1][x];
                        }
                    }
                }
            }
            k_int = k_int + newSpacing;
            kn = kn + 1;
        } while (k_int < k2);
    }

    return true;
}

bool linear2dInterpolation(dataType** originalImage, dataType** newImage, size_t imageLength, size_t imageWidth, size_t imageHeight, dataType originalSpacing, dataType newSpacing)
{
    if (originalImage == NULL || newImage == NULL)
        return false;

    size_t i, j, k, x, kn;
    dataType k_int, k1, k2;
    dataType divisionByOriginalSpacing = 1.0 / originalSpacing;

    k_int = 0; kn = 0;
    for (k = 0; k < imageHeight - 1; k++) {
        k1 = k * originalSpacing;
        k2 = k1 + originalSpacing;
        do {
            for (i = 0; i < imageLength; i++) {
                for (j = 0; j < imageWidth; j++) {
                    x = x_new(i, j, imageLength);
                    newImage[kn][x] = (dataType)(originalImage[k][x] * ((k2 - k_int) * divisionByOriginalSpacing) +
                        originalImage[k + 1][x] * ((k_int - k1) * divisionByOriginalSpacing));
                }
            }
            k_int = k_int + newSpacing;
            kn = kn + 1;
        } while (k_int < k2);
    }

    return true;
}

bool downSampling(dataType** originalImage, dataType** newImage, size_t length, size_t width, size_t height) {

    if (originalImage == NULL || newImage == NULL)
        return false;

    size_t lengthNew = (length / 2), widthNew = (width / 2), heightNew = (height / 2);

    size_t i, j, k;

    for (k = 0; k < heightNew; k++) {
        for (i = 0; i < lengthNew; i++) {
            for (j = 0; j < widthNew; j++) {
                newImage[k][x_new(i, j, lengthNew)] = originalImage[k * 2][x_new(i * 2, j * 2, length)];
            }
        }
    }
    return true;
}

bool upSampling(dataType** originalImage, dataType** newImage, size_t length, size_t width, size_t height) {

    if (originalImage == NULL || newImage == NULL)
        return false;

    size_t i, j, k, in, jn, kn, indx;

    size_t lengthNew = length * 2, widthNew = width * 2, heightNew = height * 2;

    for (k = 0; k < height; k++) {
        for (i = 0; i < length; i++) {
            for (j = 0; j < width; j++) {

                indx = x_new(i, j, length);
                in = 2 * i; jn = j * 2; kn = k * 2;

                newImage[kn][x_new(in, jn, lengthNew)] = originalImage[k][indx];
                newImage[kn][x_new(in + 1, jn, lengthNew)] = originalImage[k][indx];
                newImage[kn][x_new(in, jn + 1, lengthNew)] = originalImage[k][indx];
                newImage[kn][x_new(in + 1, jn + 1, lengthNew)] = originalImage[k][indx];

                newImage[kn + 1][x_new(in, jn, lengthNew)] = originalImage[k][indx];
                newImage[kn + 1][x_new(in + 1, jn, lengthNew)] = originalImage[k][indx];
                newImage[kn + 1][x_new(in, jn + 1, lengthNew)] = originalImage[k][indx];
                newImage[kn + 1][x_new(in + 1, jn + 1, lengthNew)] = originalImage[k][indx];
            }
        }
    }

    return true;
}

//========================================================

Point3D imageCoordToRealCoord(Point3D srcPoint, Point3D imageOrigin, VoxelSpacing imageSpacing, OrientationMatrix orientation) {
    Point3D resultPoint;
    resultPoint.x = imageOrigin.x + imageSpacing.sx * orientation.v1.x * srcPoint.x;
    resultPoint.y = imageOrigin.y + imageSpacing.sy * orientation.v2.y * srcPoint.y;
    resultPoint.z = imageOrigin.z + imageSpacing.sz * orientation.v3.z * srcPoint.z;
    return resultPoint;
}

Point3D realCoordToImageCoord(Point3D srcPoint, Point3D imageOrigin, VoxelSpacing imageSpacing, OrientationMatrix orientation) {
    Point3D resultPoint;
    resultPoint.x = (srcPoint.x - imageOrigin.x) / imageSpacing.sx * orientation.v1.x;
    resultPoint.y = (srcPoint.y - imageOrigin.y) / imageSpacing.sy * orientation.v2.y;
    resultPoint.z = (srcPoint.z - imageOrigin.z) / imageSpacing.sz * orientation.v3.z;
    return resultPoint;
}

// ======================================================

bool resizeImage(Image_Data2D oldImage, Image_Data2D newImage) {

    if (oldImage.imageDataPtr == NULL || newImage.imageDataPtr == NULL)
        return false;

    const size_t height_old = oldImage.height, width_old = oldImage.width;
    const size_t height_new = newImage.height, width_new = newImage.width;

    int i, j;
    dataType i_new, j_new;

    int i_int, j_int;

    int i_floor, i_ceil;
    int j_floor, j_ceil;

    //Compute scale factor
    dataType sx = (dataType)height_old / height_new;
    dataType sy = (dataType)width_old / width_new;

    dataType val = 0.0;

    //Map back to original coordinates
    for (i = 0; i < height_new; i++) {
        for (j = 0; j < width_new; j++) {

            //Compute corresponding point to current pixel in old image
            i_new = i * sx;
            j_new = j * sy;

            //find neighbors
            i_floor = floor(i_new);
            if (ceil(i_new) <= height_old - 1) {
                i_ceil = ceil(i_new);
            }
            else {
                i_ceil = height_old - 1;
            }

            j_floor = floor(j_new);
            if (ceil(j_new) <= width_old - 1) {
                j_ceil = ceil(j_new);
            }
            else {
                j_ceil = width_old - 1;
            }

            i_int = (int)i_new;
            j_int = (int)j_new;

            //find the closest neighbor
            if ((i_floor == i_ceil) && (j_floor == j_ceil)) {
                val = oldImage.imageDataPtr[x_new(i_int, j_int, height_old)];
            }

            if ((i_floor == i_ceil) && (j_floor != j_ceil)) {

                /*
                if (abs(i - i_ceil) > abs(i - i_floor)) {
                    val = oldImage.imageDataPtr[x_new(i_floor, j_int, height_old)];
                }
                else {
                    val = oldImage.imageDataPtr[x_new(i_ceil, j_int, height_old)];
                }
                */

                if (abs(j - j_ceil) > abs(j - j_floor)) {
                    val = oldImage.imageDataPtr[x_new(i_int, j_floor, height_old)];
                }
                else {
                    val = oldImage.imageDataPtr[x_new(i_int, j_ceil, height_old)];
                }

            }

            if (i_floor != i_ceil && (j_floor == j_ceil)) {

                /*
                if (abs(j - j_ceil) > abs(j - j_floor)) {
                    val = oldImage.imageDataPtr[x_new(i_int, j_floor, height_old)];
                }
                else {
                    val = oldImage.imageDataPtr[x_new(i_int, j_floor, height_old)];
                }
                */

                if (abs(i - i_ceil) > abs(i - i_floor)) {
                    val = oldImage.imageDataPtr[x_new(i_floor, j_int, height_old)];
                }
                else {
                    val = oldImage.imageDataPtr[x_new(i_ceil, j_int, height_old)];
                }

            }


            if (i_floor != i_ceil && (j_floor != j_ceil)) {

                dataType minDist = 100000000;

                //Top left neighbor
                dataType distTL = sqrt(pow(i - i_floor, 2) + pow(j - j_floor, 2));
                if (distTL < minDist) {
                    minDist = distTL;
                    val = oldImage.imageDataPtr[x_new(i_floor, j_floor, height_old)];
                }

                //Top right neighbor
                dataType distTR = sqrt(pow(i - i_ceil, 2) + pow(j - j_floor, 2));
                if (distTR < minDist) {
                    minDist = distTR;
                    val = oldImage.imageDataPtr[x_new(i_ceil, j_floor, height_old)];
                }

                //Bottom Left neighbor
                dataType distBL = sqrt(pow(i - i_floor, 2) + pow(j - j_ceil, 2));
                if (distBL < minDist) {
                    minDist = distBL;
                    val = oldImage.imageDataPtr[x_new(i_floor, j_ceil, height_old)];
                }

                //Bottom Right neighbor
                dataType distBR = sqrt(pow(i - i_ceil, 2) + pow(j - j_ceil, 2));
                if (distBR < minDist) {
                    minDist = distBR;
                    val = oldImage.imageDataPtr[x_new(i_ceil, j_ceil, height_old)];
                }
            }

            newImage.imageDataPtr[x_new(i, j, height_new)] = val;
        }
    }

    return true;
}

//=======================================================

Point2D getRealCoordFromImageCoord2D(Point2D srcPoint, Point2D realOrigin, PixelSpacing imageSpacing, OrientationMatrix2D orientation) {
    Point2D resultPoint;
    resultPoint.x = realOrigin.x + imageSpacing.sx * orientation.v1.x * srcPoint.x;
    resultPoint.y = realOrigin.y + imageSpacing.sy * orientation.v2.y * srcPoint.y;
    return resultPoint;
}

Point2D getImageCoordFromRealCoord2D(Point2D srcPoint, Point2D imageOrigin, PixelSpacing imageSpacing, OrientationMatrix2D orientation) {
    Point2D resultPoint;
    resultPoint.x = (srcPoint.x - imageOrigin.x) / imageSpacing.sx * orientation.v1.x;
    resultPoint.y = (srcPoint.y - imageOrigin.y) / imageSpacing.sy * orientation.v2.y;
    return resultPoint;
}

PointNeighbors2D findPointNeighbors(Point2D point, const size_t Xmax, const size_t Ymax) {

    PointNeighbors2D neighbor; 
    dataType x = point.x, y = point.y;

    int x_floor = floor(x), y_floor = floor(y), x_ceil, y_ceil;

    if (ceil(x) <= Xmax - 1) {
        x_ceil = ceil(x);
    }
    else {
        x_ceil = Xmax - 1;
    }

    if (ceil(y) <= Ymax - 1) {
        y_ceil = ceil(y);
    }
    else {
        y_ceil = Ymax - 1;
    }

    neighbor.top_left.x = x_floor;
    neighbor.top_left.y = y_floor;

    neighbor.bottom_left.x = x_floor;
    neighbor.bottom_left.y = y_ceil;

    neighbor.top_right.x = x_ceil;
    neighbor.top_right.y = y_floor;

    neighbor.bottom_right.x = x_ceil;
    neighbor.bottom_right.y = y_ceil;

    return neighbor;
}

bool ArePointsTheSame(Point2D point1, Point2D point2) {
    if ((point1.x == point2.x) && (point1.y == point2.y)) {
        return true;
    }
    else {
        return false;
    }
}

dataType getDistance2D(Point2D point1, Point2D point2) {
    return sqrt(pow(point1.x - point2.x, 2) + pow(point1.y - point2.y, 2));
}

/*
* This function return the nearest neighbor of given point
* Point2D point : given poin
*/
Point2D findNearestNeighbor2d(Point2D point, const size_t Xmax, const size_t Ymax) {

    PointNeighbors2D neighbor;
    Point2D resultPoint;

    neighbor = findPointNeighbors(point, Xmax, Ymax);

    //Case 1 :
    if (ArePointsTheSame(neighbor.top_left, neighbor.bottom_right) == true) {
        resultPoint.x = point.x;
        resultPoint.y = point.y;
    }

    //Case 2 :
    if (neighbor.top_left.x == neighbor.bottom_left.x && ArePointsTheSame(neighbor.top_left, neighbor.bottom_left) == false) {
        if (abs(point.y - neighbor.top_left.y) < abs(point.y - neighbor.bottom_left.y)) {
            resultPoint.x = point.x;
            resultPoint.y = (int)neighbor.top_left.y;
        }
        else {
            resultPoint.x = point.x;
            resultPoint.y = (int)neighbor.bottom_left.y;
        }
    }

    //Case 3 :
    if (neighbor.top_left.y == neighbor.top_right.y && ArePointsTheSame(neighbor.top_left, neighbor.top_right) == false) {
        if (abs(point.x - neighbor.top_left.x) < abs(point.x - neighbor.top_right.x)) {
            resultPoint.x = (int)neighbor.top_left.x;
            resultPoint.y = point.y;
        }
        else {
            resultPoint.x = (int)neighbor.top_left.x;
            resultPoint.y = point.y;
        }
    }

    //Case 4 :
    if (ArePointsTheSame(neighbor.top_left, neighbor.bottom_right) == false) {

        dataType distance[4];
        distance[0] = getDistance2D(point, neighbor.top_left);
        distance[1] = getDistance2D(point, neighbor.top_right);
        distance[2] = getDistance2D(point, neighbor.bottom_left);
        distance[3] = getDistance2D(point, neighbor.bottom_right);

        dataType min_dist = 1000000000000;
        int k, min_indice;

        for (k = 0; k < 4; k++) {
            if (distance[k] < min_dist) {
                min_dist = distance[k];
                min_indice = k;
            }
        }

        if (min_indice == 0) {
            resultPoint.x = (int)neighbor.top_left.x;
            resultPoint.y = (int)neighbor.top_left.y;
        }
        if (min_indice == 1) {
            resultPoint.x = (int)neighbor.top_right.x;
            resultPoint.y = (int)neighbor.top_right.y;
        }
        if (min_indice == 2) {
            resultPoint.x = (int)neighbor.bottom_left.x;
            resultPoint.y = (int)neighbor.bottom_left.y;
        }
        if (min_indice == 3) {
            resultPoint.x = (int)neighbor.bottom_right.x;
            resultPoint.y = (int)neighbor.bottom_right.y;
        }

    }

    return resultPoint;
}

bool resizeImageFromImageCoordToRealCoord(Image_Data2D src_image, Image_Data2D dest_image, OrientationMatrix2D orientation) {
    
    if (src_image.imageDataPtr == NULL || dest_image.imageDataPtr == NULL)
        return false;

    const size_t src_height = src_image.height, src_width = src_image.width;
    const size_t dest_height = dest_image.height, dest_width = dest_image.width;

    int i, j, i_int, j_int;
    Point2D src_point, dest_point, int_point;
    PointNeighbors2D neighbor_points;

    //indexes for neighbors 
    int i_floor, i_ceil, j_floor, j_ceil;

    dataType pixel_value = 0.0;

    for (i = 0; i < dest_height; i++) {
        for (j = 0; j < dest_width; j++) {

            dest_point.x = i;
            dest_point.y = j;

            //find the corresponding point to the real CS points 
            src_point = getImageCoordFromRealCoord2D(dest_point, dest_image.origin, dest_image.spacing, orientation);

            //collection of indices
            //please separate to specific function 

            /*
            i_floor = floor(src_point.x);
            if (ceil(src_point.x) <= src_height - 1) {
                i_ceil = ceil(src_point.x);
            }
            else {
                i_ceil = src_height - 1;
            }

            j_floor = floor(src_point.y);
            if (ceil(src_point.y) <= src_width - 1) {
                j_ceil = ceil(src_point.y);
            }
            else {
                j_ceil = src_width - 1;
            }
            */

            //Find current point neighbors
            //neighbor_points = findPointNeighbors(src_point, src_height, src_width);

            //i_int = (int)src_point.x;
            //j_int = (int)src_point.y;

            //specific interpolation method
            // please extract it to specific function
            //find the closest neighbor

            //if ((i_floor == i_ceil) && (j_floor == j_ceil)) {
            //    pixel_value = src_image.imageDataPtr[x_new(i_int, j_int, src_height)];
            //}

            /*
            if (ArePointsTheSame(neighbor_points.top_left, neighbor_points.bottom_right) == true) {
                pixel_value = src_image.imageDataPtr[x_new(i_int, j_int, src_height)];
            }
            */

            /*
            if ((i_floor == i_ceil) && (j_floor != j_ceil)) {

                if (abs(i - i_ceil) > abs(i - i_floor)) {
                    pixel_value = src_image.imageDataPtr[x_new(i_floor, j_int, src_height)];
                }
                else {
                    pixel_value = src_image.imageDataPtr[x_new(i_ceil, j_int, src_height)];
                }
            }
            */

            /*
            if (neighbor_points.top_left.x == neighbor_points.bottom_left.x && ArePointsTheSame(neighbor_points.top_left, neighbor_points.bottom_left) == false) {
                if (abs(j - neighbor_points.top_left.y) < abs(j - neighbor_points.bottom_left.y)) {
                    pixel_value = pixel_value = src_image.imageDataPtr[x_new(i_int, (int)neighbor_points.top_left.y, src_height)];
                }
                else {
                    pixel_value = pixel_value = src_image.imageDataPtr[x_new(i_int, (int)neighbor_points.bottom_left.y, src_height)];
                }
            }
            */

            /*
            if (i_floor != i_ceil && (j_floor == j_ceil)) {

                if (abs(j - j_ceil) > abs(j - j_floor)) {
                    pixel_value = src_image.imageDataPtr[x_new(i_int, j_floor, src_height)];
                }
                else {
                    pixel_value = src_image.imageDataPtr[x_new(i_int, j_floor, src_height)];
                }
            }
            */

            /*
            if (neighbor_points.top_left.y == neighbor_points.top_right.y && ArePointsTheSame(neighbor_points.top_left, neighbor_points.top_right) == false) {
                if (abs(i - neighbor_points.top_left.x) < abs(i - neighbor_points.top_right.x)) {
                    pixel_value = src_image.imageDataPtr[x_new((int)neighbor_points.top_left.x, j_int, src_height)];
                }
                else {
                    pixel_value = src_image.imageDataPtr[x_new((int)neighbor_points.top_right.x, j_int, src_height)];
                }
            }
            */

            /*
            if (i_floor != i_ceil && (j_floor != j_ceil)) {

                dataType minDist = 100000000;

                //Top left neighbor
                dataType distTL = sqrt(pow(i - i_floor, 2) + pow(j - j_floor, 2));
                if (distTL < minDist) {
                    minDist = distTL;
                    pixel_value = src_image.imageDataPtr[x_new(i_floor, j_floor, src_height)];
                }

                //Top right neighbor
                dataType distTR = sqrt(pow(i - i_ceil, 2) + pow(j - j_floor, 2));
                if (distTR < minDist) {
                    minDist = distTR;
                    pixel_value = src_image.imageDataPtr[x_new(i_ceil, j_floor, src_height)];
                }

                //Bottom Left neighbor
                dataType distBL = sqrt(pow(i - i_floor, 2) + pow(j - j_ceil, 2));
                if (distBL < minDist) {
                    minDist = distBL;
                    pixel_value = src_image.imageDataPtr[x_new(i_floor, j_ceil, src_height)];
                }

                //Bottom Right neighbor
                dataType distBR = sqrt(pow(i - i_ceil, 2) + pow(j - j_ceil, 2));
                if (distBR < minDist) {
                    minDist = distBR;
                    pixel_value = src_image.imageDataPtr[x_new(i_ceil, j_ceil, src_height)];
                }
            }
            */

            /*
            if (ArePointsTheSame(neighbor_points.top_left, neighbor_points.bottom_right) == false) {

                dataType distance[4];
                distance[0] = getDistance2D(dest_point, neighbor_points.top_left);
                distance[1] = getDistance2D(dest_point, neighbor_points.top_right);
                distance[2] = getDistance2D(dest_point, neighbor_points.bottom_left);
                distance[3] = getDistance2D(dest_point, neighbor_points.bottom_right);

                dataType min_dist = 1000000000000;
                int k, min_indice;

                for (k = 0; k < 4; k++) {
                    if (distance[k] < min_dist) {
                        min_dist = distance[k];
                        min_indice = k;
                    }
                }

                if (min_indice == 0) {
                    pixel_value = src_image.imageDataPtr[x_new((int)neighbor_points.top_left.x, (int)neighbor_points.top_left.y, src_height)];
                }
                if (min_indice == 1) {
                    pixel_value = src_image.imageDataPtr[x_new((int)neighbor_points.top_right.x, (int)neighbor_points.top_right.y, src_height)];
                }
                if (min_indice == 2) {
                    pixel_value = src_image.imageDataPtr[x_new((int)neighbor_points.bottom_left.x, (int)neighbor_points.bottom_left.y, src_height)];
                }
                if (min_indice == 3) {
                    pixel_value = src_image.imageDataPtr[x_new((int)neighbor_points.bottom_right.x, (int)neighbor_points.bottom_right.y, src_height)];
                }

            }
            */

            int_point = findNearestNeighbor2d(src_point, src_height, src_width);

            //dest_image.imageDataPtr[x_new(i, j, dest_height)] = pixel_value;

            dest_image.imageDataPtr[x_new(i, j, dest_height)] = src_image.imageDataPtr[x_new((int)int_point.x, (int)int_point.y, src_height)];

        }
    }
    return true;
}
