#if !defined(__EXTRACT_H__)
#define __EXTRACT_H__

class RCFile_C;

int checkBitmaps(RCFile_C *pRC, const wchar_t *pszBaseDir);
int checkCursors(RCFile_C *pRC, const wchar_t *pszBaseDir);

#endif //!defined(__EXTRACT_H__)
