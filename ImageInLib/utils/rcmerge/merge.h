#if !defined(__MERGE_H__)
#define __MERGE_H__

class RCFile_C;

int mergeTexts(RCFile_C *pSourceRC, RCFile_C *pTranslatedRC, const wchar_t *pszSourceRCFile, const wchar_t *pszOutputRCFile);

#endif //!defined(__MERGE_H__)
