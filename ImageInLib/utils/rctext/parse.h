#if !defined(__PARSE_H__)
#define __PARSE_H__

class RCFile_C;
class Resource_C;
class MenuPopup_C;
class VersionInfoBlock_C;
class StringTable_C;

typedef struct yy_parse_state {

    const wchar_t *fileName;
	wchar_t *scanBuffer;
	size_t scanBufferLength;
	RCFile_C   *rc;

	// parsed file position
    int curLine;
    int curColumn;

	// parser state
	StringTable_C *stringtable;
	Resource_C *resource;
	std::wstring resourceName;
	std::wstring resourceType;
	std::wstring resourceOptions;
	std::wstring resourceItem;
	std::wstring controlID;
	std::wstring controlText;
	std::wstring controlType;
	std::wstring controlOptions;
	std::wstring preprocessorParams;
	MenuPopup_C *popup;
	VersionInfoBlock_C *infoblock;
	int conditionState;
	stack<int> conditionStack;

	// scanner state
    int ctxComment;
    int ctxDirective;
    int trigSkip;

} yy_parse_state;


void parseError(yy_parse_state *yystate, const char *pszErrorMsg);
int  checkCondition(const wchar_t *pszCondition, bool bComplex);

RCFile_C *parseRCFile(const wchar_t *pszFileName, bool bRequireUnicode);


#endif //!defined(__PARSE_H__)
