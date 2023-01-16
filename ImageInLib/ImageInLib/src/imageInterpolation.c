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
    return false;
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