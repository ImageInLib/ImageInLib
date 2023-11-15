#include <memory.h>
#include <corecrt_malloc.h>

#include "solvers.h"

bool sherman_morris(dataType* a, dataType* b, dataType* c, dataType alpha, dataType beta, dataType* x, dataType* ps, const size_t N)
{
    if (a == NULL || b == NULL || c == NULL || x == NULL || ps == NULL) {
        return false;
    }

    int i = 0;
    dataType fact = 0, gamma = 0;
    dataType* bb = (dataType*)malloc(sizeof(dataType*) * (N + 2));
    if (bb == NULL)
    {
        return false;
    }
    dataType* u = (dataType*)malloc(sizeof(dataType*) * (N + 2));
    if (u == NULL) {
        free(bb);
        return false;
    }
    dataType* z = (dataType*)malloc(sizeof(dataType*) * (N + 2));
    if (z == NULL) {
        free(bb);
        free(u);
        return false;
    }

    memcpy(bb, b, sizeof(dataType) * (N + 2));
    gamma = -b[1];
    bb[1] -= gamma;
    bb[N] -= alpha * beta / gamma;

    thomas(a, bb, c, x, ps, N);

    memset(u, 0, sizeof(dataType) * (N + 2));

    u[1] = gamma;
    u[N] = alpha;

    thomas(a, bb, c, z, u, N);

    fact = (x[1] + beta * x[N] / gamma) / ((dataType)1.0 + z[1] + beta * z[N] / gamma);

    for (i = 1; i < N + 1; i++)
    {
        x[i] -= fact * z[i];
    }

    free(bb);
    free(u);
    free(z);

    return true;
}

bool thomas(dataType* a, dataType* b, dataType* c, dataType* x, dataType* ps, const size_t N)
{
    if (a == NULL || b ==NULL || c == NULL || x == NULL || ps == NULL) {
        return false;
    }

    size_t i = 0;
    dataType m = 0;
    dataType* bb = (dataType*)malloc(sizeof(dataType*) * (N + 2));
    if (bb == NULL) {
        return false;
    }
    dataType* pps = (dataType*)malloc(sizeof(dataType*) * (N + 2));
    if (pps == NULL) {
        free(bb);
        return false;
    }

    memcpy(bb, b, sizeof(dataType) * (N + 2));
    memcpy(pps, ps, sizeof(dataType) * (N + 2));

    //Forward elimination phase

    for (i = 2; i < N + 1; i++)
    {
        m = a[i] / bb[i - 1];
        bb[i] -= m * c[i - 1];
        pps[i] -= m * pps[i - 1];
    }

    //Bakward substitution phase

    x[N] = pps[N] / bb[N];

    for (i = N - 1; i >= 1; i--)
    {
        x[i] = (pps[i] - c[i] * x[i + 1]) / bb[i];
    }

    free(bb);
    free(pps);

    return true;
}
