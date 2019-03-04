#pragma once
#include "common_functions.h"
#include "registration_params.h"

void imageRegistration(dataType **, dataType **, dataType **, const size_t imageHeight, const size_t imageLength, const size_t imageWidth, registrationParams params, optimizationMethod optimization);