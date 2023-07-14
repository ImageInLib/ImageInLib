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

/*
 * Get 3D point in real world from 3D image point
 * Three operations are done here : image scaling, image translation and image rotation
*/
Point3D imageCoordToRealCoord(Point3D srcPoint, Point3D imageOrigin, VoxelSpacing imageSpacing, OrientationMatrix orientation) {
    Point3D resultPoint;
    resultPoint.x = imageOrigin.x + imageSpacing.sx * (orientation.v1.x + orientation.v1.y + orientation.v1.z) * srcPoint.x;
    resultPoint.y = imageOrigin.y + imageSpacing.sy * (orientation.v2.x + orientation.v2.y + orientation.v2.z) * srcPoint.y;
    resultPoint.z = imageOrigin.z + imageSpacing.sz * (orientation.v3.x + orientation.v3.y + orientation.v3.z) * srcPoint.z;
    return resultPoint;
}


//Get image coordinate from real coordinate
//TO DO : find manually the matrix inverse in order to write the final function 
Point3D realCoordToImageCoord(Point3D srcPoint, Point3D realOrigin, VoxelSpacing imageSpacing, OrientationMatrix orientation) {
    Point3D resultPoint;

    return resultPoint;
}

/*
* This function perform 2d nearest neighbor interpolation (expansion or shrinking)
*/
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

            //Compute correponding point to current pixel in old image
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

                if (abs(i - i_ceil) > abs(i - i_floor)) {
                    val = oldImage.imageDataPtr[x_new(i_floor, j_int, height_old)];
                }
                else {
                    val = oldImage.imageDataPtr[x_new(i_ceil, j_int, height_old)];
                }

            }
            if (i_floor != i_ceil && (j_floor == j_ceil)) {

                if (abs(j - j_ceil) > abs(j - j_floor)) {
                    val = oldImage.imageDataPtr[x_new(i_int, j_floor, height_old)];
                }
                else {
                    val = oldImage.imageDataPtr[x_new(i_int, j_floor, height_old)];
                }
            }
            if (i_floor != i_ceil && (j_floor != j_ceil)) {

                dataType minDist = 100000000;

                //Top left neigh
                dataType distTL = sqrt(pow(i - i_floor, 2) + pow(j - j_floor, 2));
                if (distTL < minDist) {
                    minDist = distTL;
                    val = oldImage.imageDataPtr[x_new(i_floor, j_floor, height_old)];
                }

                //Top right neigh
                dataType distTR = sqrt(pow(i - i_ceil, 2) + pow(j - j_floor, 2));
                if (distTR < minDist) {
                    minDist = distTR;
                    val = oldImage.imageDataPtr[x_new(i_ceil, j_floor, height_old)];
                }

                //Bottom Left neigh
                dataType distBL = sqrt(pow(i - i_floor, 2) + pow(j - j_ceil, 2));
                if (distBL < minDist) {
                    minDist = distBL;
                    val = oldImage.imageDataPtr[x_new(i_floor, j_ceil, height_old)];
                }

                //Bottom Right neigh
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
