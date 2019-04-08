#pragma once

#include <atlstr.h>

namespace DataModel {

	class CFunction
	{
		public:
			CFunction(UINT nThreadID)
			{
				Start = -1;
				Stop = -1;
				Level = 0;
				ThreadID = nThreadID;
			}
		// members
		public:
			double  Start;
			double  Stop;
			CString Description;
			UINT    Level;
			UINT    ThreadID;
	};
}
