#include "stdafx.h"
#include "FormatTime.h"
#include <algorithm>

CString FormatTimeInterval(double dTime, int precision)
{
	CString str;
	if (dTime >= 1)    // display seconds
	{
		str.Format(_T("%0.*f s"), precision, dTime );
	}
	else if (dTime >= 1e-3 || dTime == 0.0)
	{
		str.Format(_T("%0.*f ms"), precision, dTime * 1e3);
	}
	else if (dTime >= 1e-6)    // display microseconds
	{
		str.Format(_T("%0.*f 탎"), precision, dTime * 1e6);
	}
	else    // display nanoseconds
	{
		str.Format(_T("%0.*f ns"), precision, dTime * 1e9);
	}
	return str;
}

CString FormatAbsoluteTime(double dTime, int precisionAfterSeconds)
{
	CString str;
	if (dTime >= 1)    // display seconds
	{
		str.Format(_T("%0.*f s"), precisionAfterSeconds ? precisionAfterSeconds : 1, dTime);
	}
	else if (dTime >= 1e-3 || dTime == 0.0)
	{
		precisionAfterSeconds = std::max<int>(precisionAfterSeconds - 3, 1);
		str.Format(_T("%0.*f ms"), precisionAfterSeconds ? precisionAfterSeconds : 1, dTime * 1e3);
	}
	else if (dTime >= 1e-4)    // display microseconds
	{
		precisionAfterSeconds = std::max<int>(precisionAfterSeconds - 6, 1);
		str.Format(_T("%0.*f 탎"), precisionAfterSeconds ? precisionAfterSeconds : 1, dTime * 1e6);
	}
	else if (dTime >= 1e-5)
	{
		precisionAfterSeconds = std::max<int>(precisionAfterSeconds - 6, 1);
		str.Format(_T("%0.*f 탎"), precisionAfterSeconds ? precisionAfterSeconds : 2, dTime * 1e6);
	}
	else if (dTime >= 1e-6)
	{
		precisionAfterSeconds = std::max<int>(precisionAfterSeconds - 6, 1);
		str.Format(_T("%0.*f 탎"), precisionAfterSeconds ? precisionAfterSeconds : 3, dTime * 1e6);
	}
	else    // display nanoseconds
	{
		precisionAfterSeconds = std::max<int>(precisionAfterSeconds - 9, 1);
		str.Format(_T("%0.*f ns"), precisionAfterSeconds ? precisionAfterSeconds : 1, dTime * 1e9);
	}
	return str;
}
