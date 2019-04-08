// Security.cpp
//
//////////////////////////////////////////////////////////////////////

#include "headers.h"
#include "Security.h"

using namespace ::Security;


bool Security::operator<(const CBinary &l, const CBinary &r)
{
	int c = memcmp((PBYTE)l, (PBYTE)r, min(l.size(), r.size())); 
	if (c < 0 || (c == 0 && l.size() < r.size())) 
		return true; 
	return false; 
}

bool Security::operator==(const CBinary &l, const CBinary &r)
{ 
	if (l.size() != r.size())
		return false;
	if (memcmp((PBYTE)l, (PBYTE)r, l.size()) == 0)
		return true;
	return false; 
}

bool Security::operator!=(const CBinary &l, const CBinary &r)
{ 
	if (l.size() != r.size())
		return true;
	if (memcmp((PBYTE)l, (PBYTE)r, l.size()) == 0)
		return false;
	return true; 
}
