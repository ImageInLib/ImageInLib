#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include "common_functions.h"
#include "../src/registration_params.h"

	void imageRegistration(dataType **, dataType **, dataType **, const size_t imageHeight, const size_t imageLength, const size_t imageWidth, Registration_Params params, Optimization_Method optimization);

#ifdef __cplusplus
}
#endif