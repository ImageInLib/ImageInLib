#ifndef __ReducedColorsCursor_h__INCLUDED__
#define __ReducedColorsCursor_h__INCLUDED__

#include "mutex.h"
#include <map>

class ReducedColorsCursor
{
public:
	static HCURSOR LoadCursor(UINT nIDResource, int bits);
	
protected:
	static mutex m_cursorsLock;

	struct CURSOR_INFO
	{
		HCURSOR m_hCursor;
		HBITMAP m_hbmMask;
		HBITMAP m_hbmColor;
	};	
	static std::map<UINT, CURSOR_INFO> m_cursors1BPP;
	static std::map<UINT, CURSOR_INFO> m_cursors4BPP;
	static std::map<UINT, CURSOR_INFO> m_cursors8BPP;
	static std::map<UINT, CURSOR_INFO> m_cursors24BPP;
};

#endif //__ReducedColorsCursor_h__INCLUDED__
