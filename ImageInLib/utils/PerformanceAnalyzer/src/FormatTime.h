#pragma once

#include <atlstr.h>

static const int DefaultIntervalPrecision = 3;
static const int DefaultPrecisionAfterSeconds = 7;

CString FormatTimeInterval(double dTime, int precision);
CString FormatAbsoluteTime(double dTime, int precisionAfterSeconds = 0);