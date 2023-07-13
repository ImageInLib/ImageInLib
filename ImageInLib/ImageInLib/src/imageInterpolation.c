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
                        if ( (k_int - k1) < (k2 - k_int)) {
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
 * Get real world cordinate from image coordinate
 * srcPoint : contains the voxel indexes
 * realOrigin : image origin in real world
 * imageSpacing : - distance between voxels in x and y direction
 *                - distance between slices
 * orientation : IJK(1,1,1), RAS(-1,-1,1) or LPS(1,-1,-1)
*/
Point3D imageCoordToRealCoord(imgPoint3D srcPoint, Point3D realOrigin, Spacing3D imageSpacing, Point3D orientation) { 
    Point3D resultPoint;
    resultPoint.x = realOrigin.x + orientation.x * imageSpacing.sx * srcPoint.x;
    resultPoint.y = realOrigin.y + orientation.y * imageSpacing.sy * srcPoint.y;
    resultPoint.z = realOrigin.z + orientation.z * imageSpacing.sz * srcPoint.z;
    return resultPoint;
}

//Get image coordinate from real coordinate
imgPoint3D realCoordToImageCoord(Point3D srcPoint, Point3D realOrigin, Spacing3D imageSpacing, Point3D orientation) {
    imgPoint3D resultPoint;
    resultPoint.x = (dataType)((srcPoint.x - realOrigin.x) / (orientation.x * imageSpacing.sx));
    resultPoint.y = (dataType)((srcPoint.y - realOrigin.y) / (orientation.y * imageSpacing.sy));
    resultPoint.z = (dataType)((srcPoint.z - realOrigin.z) / (orientation.z * imageSpacing.sz));
    return resultPoint;
}

bool resizeImage(dataType* oldImage, dataType* newImage) {

    if (oldImage == NULL || newImage == NULL)
        return false;

    const size_t height_old = 512, width_old = 512;
    const size_t height_new = 600, width_new = 600;

    size_t i, j;
    dataType i_new, j_new;

    size_t i_int, j_int;

    size_t i_floor, i_ceil;
    size_t j_floor, j_ceil;

    //Compute scale factor
    dataType sx = height_old / height_new;
    dataType sy = width_old / width_new;

    dataType val = 0.0;

    //Map back to original coordinates
    for (i = 0; i < height_new; i++) {
        for (j = 0; j < width_new; j++) {

            i_new = i * sx;
            j_new = j * sy;

            i_floor = floor(i_new); 
            if (ceil(i_new) < height_old - 1) {
                i_ceil = ceil(i_new);
            }
            else {
                i_ceil = height_old - 1;
            }

            j_floor = floor(j_new);
            if (ceil(j_new) < width_old - 1) {
                j_ceil = ceil(j_new);
            }
            else {
                j_ceil = width_old - 1;
            }

            i_int = (size_t)i_new;
            j_int = (size_t)j_new;

            if (i_floor == i_ceil && j_floor == j_ceil) {
                val = oldImage[x_new(i_int, j_int, height_old)];
            }
            else {
                if (i_floor == i_ceil) {
                    dataType val1 = oldImage[x_new(i_int, j_floor, height_old)];
                    dataType val2 = oldImage[x_new(i_int, j_ceil, height_old)];
                    val = val1 * (j_ceil - j) + val2 * (j - j_floor);
                }
                else {
                    if (j_floor == j_ceil) {
                        dataType val1 = oldImage[x_new(i_floor, j_int, height_old)];
                        dataType val2 = oldImage[x_new(i_ceil, j_int, height_old)];
                        val = val1 * (i_ceil - i) + val2 * (i - i_floor);
                    }
                    else {
                        dataType val1 = oldImage[x_new(i_floor, j_floor, height_old)];
                        dataType val2 = oldImage[x_new(i_ceil, j_floor, height_old)];
                        dataType val3 = oldImage[x_new(i_floor, j_ceil, height_old)];
                        dataType val4 = oldImage[x_new(i_ceil, j_ceil, height_old)];
                        
                        dataType val11 = val1 * (i_ceil - i) + val2 * (i - i_floor);
                        dataType val12 = val3 * (i_ceil - i) + val4 * (i - i_floor);
                        
                        val = val11 * (j_ceil - j) + val12 * (j - j_floor);
                    }
                }
            }
            newImage[x_new(i, j, height_new)] = val;
        }
    }

    return true;
}
