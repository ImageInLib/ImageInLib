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
// !!!!!!!!!! Not done : To be fixed !!!!!!!!!!!!!!!!!
bool interpolateToRealDimension(patientImageData imageSrc, const char* outputPathPtr) {
    
    if (imageSrc.dataPtr == NULL)
        return false;

    size_t k, i, j;
    dataType k_real = 0, i_real = 0, j_real = 0;
    size_t k_int = 0, i_int = 0, j_int = 0;

    size_t height = imageSrc.height, length = imageSrc.length, width = imageSrc.width;
    const size_t real_height = (size_t)(height * imageSrc.toRealCoordinates.sz);
    const size_t real_length = (size_t)(length * imageSrc.toRealCoordinates.sx);
    const size_t real_width = (size_t)(width * imageSrc.toRealCoordinates.sy);

    dataType** destImage = (dataType**)malloc(sizeof(dataType*) * real_height);
    for (k = 0; k < real_height; k++) {
        destImage[k] = (dataType*)malloc(sizeof(dataType) * real_length * real_width);
    }
    if (destImage == NULL)
        return false;

    for (k = 0; k < height; k++) {

        k_real = k * imageSrc.toRealCoordinates.sz;
        dataType diffZ1 = abs(k_real - k);
        dataType diffZ2 = abs(k_real - (k + 1));
        if (diffZ1 < diffZ2) {
            k_int = k;
        }
        else {
            k_int = k + 1;
        }

        for (i = 0; i < length; i++) {

            i_real = i * imageSrc.toRealCoordinates.sx;
            dataType diffX1 = abs(i_real - i);
            dataType diffX2 = abs(i_real - (i + 1));
            if (diffX1 < diffX2) {
                i_int = i;
            }
            else {
                i_int = i + 1;
            }

            for (j = 0; j < width; j++) {

                j_real = j * imageSrc.toRealCoordinates.sy;
                dataType diffY1 = abs(j_real - j);
                dataType diffY2 = abs(j_real - (j + 1));
                if (diffY1 < diffY2) {
                    j_int = j;
                }
                else {
                    j_int = j + 1;
                }
                destImage[k_int][x_new(i_int, j_int, real_length)] = imageSrc.dataPtr[k][x_new(i, j, length)];
            }
        }
    }

    printf("The dimensions of new image are \n Length = %zd, Width = %zd, Heigth = %zd", real_length, real_width, real_height);
    Storage_Flags flags = { false, false };
    store3dDataArrayD(destImage, real_length, real_width, real_height, outputPathPtr, flags);
    
    for (k = 0; k < real_height; k++) {
        free(destImage[k]);
    }
    free(destImage);

    return true;
}
*/

bool imageCoordToRealCoord(Point3D srcPoint, Spacing imageSpacing, Point3D destPoint) {
    destPoint.x = srcPoint.x * imageSpacing.sx;
    destPoint.y = srcPoint.y * imageSpacing.sy;
    destPoint.z = srcPoint.z * imageSpacing.sz;
}

bool IJKtoRAS(Point3D srcPoint, Spacing imageSpacing, Point3D destPoint) {

    destPoint.x = -srcPoint.x + imageSpacing.sx;
    destPoint.y = -srcPoint.y + imageSpacing.sy;
    destPoint.z = srcPoint.z + imageSpacing.sz;

    return true;
}

bool IJKtoLPS(Point3D srcPoint, Spacing imageSpacing, Point3D destPoint) {

    destPoint.x = srcPoint.x + imageSpacing.sx;
    destPoint.y = -srcPoint.y + imageSpacing.sy;
    destPoint.z = -srcPoint.z + imageSpacing.sz;

    return true;
}

//2D images for test

bool interpolateToRealDimension2D(patientImageData2D imageSrc, Spacing2D newSpacing, const char* outputPathPtr) {

    size_t i, j;

    size_t height = imageSrc.height, width = imageSrc.width;
    size_t height_new = (size_t)(height * imageSrc.toRealCoordinates.sx / newSpacing.sx);
    size_t width_new = (size_t)(width * imageSrc.toRealCoordinates.sy / newSpacing.sy);

    dataType* destImage = (dataType*)malloc(sizeof(dataType) * height_new * width_new);
    if (destImage == NULL)
        return false;

    //Interpolation in x direction
    dataType i1, i2, i_int = 0;
    size_t i_new = 0;
    for (i = 0; i < height - 1; i++) {
        i1 = imageSrc.toRealCoordinates.sx * i;
        i2 = imageSrc.toRealCoordinates.sx * (i + 1);
        do {
            for (j = 0; j < width; j++) {
                if (i_int == 0) {
                    destImage[x_new(i_new, j, height_new)] = imageSrc.dataPtr[x_new(i, j, height)];
                }
                else {
                    if ((i_int - i1) < (i2 - i_int)) {
                        destImage[x_new(i_new, j, height_new)] = imageSrc.dataPtr[x_new(i, j, height)];
                    }
                    else {
                        destImage[x_new(i_new, j, height_new)] = imageSrc.dataPtr[x_new(i + 1, j, height)];
                    }
                }
            }
            i_int = i_int + newSpacing.sx;
            i_new = i_new + 1;
        } while (i_int < i);
    }

    //Interpolation in y direction
    dataType j1, j2, j_int = 0;
    size_t j_new = 0;
    for (j = 0; j < width - 1; j++) {
        j1 = imageSrc.toRealCoordinates.sy * j;
        j2 = imageSrc.toRealCoordinates.sy * (j + 1);
        do {
            for (i = 0; i < height; i++) {
                if (j_int == 0) {
                    destImage[x_new(i, j_new, height_new)] = imageSrc.dataPtr[x_new(i, j, height)];
                }
                else {
                    if ((j_int - j1) < (j2 - j_int)) {
                        destImage[x_new(i, j_new, height_new)] = imageSrc.dataPtr[x_new(i, j, height)];
                    }
                    else {
                        destImage[x_new(i, j_new, height_new)] = imageSrc.dataPtr[x_new(i, j + 1, height)];
                    }
                }
            }
            j_int = j_int + newSpacing.sy;
            j_new = j_new + 1;
        } while (j_int < j);
    }

    printf("The dimensions of new image are \n Heigth = %zd, Width = %zd", height_new, width_new);
    Storage_Flags flags = { false, false };
    store2dRawData(destImage, height_new, width_new, outputPathPtr, flags);

    free(destImage);

    return true;
}
