%{
// FLEX lexical analyzer for parsing .RC files
//
// flex -B -U -oscanner.cpp scanner.lex
//
// Tested with flex version 2.5.4
//

#include "headers.h"
#include "grammar.h"
#ifndef fileno
#define fileno _fileno
#endif
#ifndef read
#define read _read
#endif
#ifndef isatty
#define isatty _isatty
#endif
#include "parse.h"

#define YY_INPUT(buf,result,max_size) result = 0;

static void count(yy_parse_state *yyinput);
int yylex(int *yylval, void *yyinput);

#define YY_DECL int _yylex(int *yylval, void *yyinput)
#define yystate ((yy_parse_state*)yyinput)
#define yyerror(msg) parseError((yy_parse_state*)yyinput, (msg))

#define COUNT count(yystate)

%}

%x _comment_
%x _param_
%x _skip_

D			[0-9]
L			[a-zA-Z_]
H			[0-9a-fA-F]
IS			(u|U|l|L)*

%%
	if (yystate->trigSkip == 1) {
		yystate->trigSkip = 0;
		BEGIN(_skip_);
	} else if (yystate->trigSkip == -1) {	
		yystate->trigSkip = 0;
		BEGIN(INITIAL);
	}

	/*
	  =========================
              comments
      =========================
	*/

"/*"			{ COUNT; yystate->ctxComment = YY_START; BEGIN(_comment_); }
<_comment_>{
[^*\n]*
[^*\n]*\n	   	
"*"+[^*/\n]*
"*"+[^*/\n]*\n	
"*"+"/"	   		{ COUNT; BEGIN(yystate->ctxComment); }
}
"//".*\n		{ COUNT; }
";".*\n			{ COUNT; }

	/*
	  =========================
	   preprocessor directives
	  =========================
	*/

<_skip_>.               { COUNT; }
<_skip_>[ \t\v\n\f]     { COUNT; }

<INITIAL,_skip_>{
^[ \t]*"#define"[ \t]+	{ COUNT; yystate->ctxDirective = YY_START; BEGIN(_param_); return(DEFINE_DIRECTIVE); }
^[ \t]*"#elif"[ \t]+	{ COUNT; yystate->ctxDirective = YY_START; BEGIN(_param_); return(ELIF_DIRECTIVE); }
^[ \t]*"#else"[ \t]*	{ COUNT; return(ELSE_DIRECTIVE); }
^[ \t]*"#endif"[ \t]*	{ COUNT; return(ENDIF_DIRECTIVE); }
^[ \t]*"#if"[ \t]+		{ COUNT; yystate->ctxDirective = YY_START; BEGIN(_param_); return(IF_DIRECTIVE); }
^[ \t]*"#ifdef"[ \t]+	{ COUNT; yystate->ctxDirective = YY_START; BEGIN(_param_); return(IFDEF_DIRECTIVE); }
^[ \t]*"#ifndef"[ \t]+	{ COUNT; yystate->ctxDirective = YY_START; BEGIN(_param_); return(IFNDEF_DIRECTIVE); }
^[ \t]*"#include"[ \t]+	{ COUNT; yystate->ctxDirective = YY_START; BEGIN(_param_); return(INCLUDE_DIRECTIVE); }
^[ \t]*"#undef"[ \t]+	{ COUNT; yystate->ctxDirective = YY_START; BEGIN(_param_); return(UNDEF_DIRECTIVE); }
^[ \t]*"#pragma"[ \t]+	{ COUNT; yystate->ctxDirective = YY_START; BEGIN(_param_); return(PRAGMA_DIRECTIVE); }
}
<_param_>[ \t]          { COUNT; }
<_param_>"/*"			{ COUNT; yystate->ctxComment = YY_START; BEGIN(_comment_); }
<_param_>\\\n			{ COUNT; }
<_param_>\n				{ COUNT; BEGIN(yystate->ctxDirective); return(DIRECTIVE_END); }
<_param_>[^ "\t\\/\n]*	{ COUNT; return(DIRECTIVE_PARAM); }
<_param_>\"(\\.|[^\\"])*\" { COUNT; return(DIRECTIVE_PARAM); }
<_param_>"//".*[^\n]	{ COUNT; }

	/*
	  =========================
	          keywords
	  =========================
	*/
("BEGIN"|"{")		{ COUNT; return(BLOCK_BEGIN); }
("END"|"}")			{ COUNT; return(BLOCK_END); }

"ACCELERATORS"		{ COUNT; return(STM_ACCELERATORS); }
"BITMAP"		    { COUNT; return(STM_BITMAP); }
"CURSOR"		    { COUNT; return(STM_CURSOR); }
"DIALOG"		    { COUNT; return(STM_DIALOG); }
"DIALOGEX"		    { COUNT; return(STM_DIALOGEX); }
"FONT"		    	{ COUNT; return(STM_FONT); }
"ICON"			    { COUNT; return(STM_ICON); }
"MENU"			    { COUNT; return(STM_MENU); }
"MENUEX"		    { COUNT; return(STM_MENUEX); }
"MESSAGETABLE"		{ COUNT; return(STM_MESSAGETABLE); }
"RCDATA"			{ COUNT; return(STM_RCDATA); }
"STRINGTABLE"		{ COUNT; return(STM_STRINGTABLE); }
"VERSIONINFO"		{ COUNT; return(STM_VERSIONINFO); }
"LANGUAGE"			{ COUNT; return(STM_LANGUAGE); }

"DESIGNINFO"		{ COUNT; return(DESIGNINFO); }
"TOOLBAR"			{ COUNT; return(STM_TOOLBAR); }

"AUTO3STATE"		{ COUNT; return(AUTO3STATE); }
"AUTOCHECKBOX"		{ COUNT; return(AUTOCHECKBOX); }
"AUTORADIOBUTTON"	{ COUNT; return(AUTORADIOBUTTON); }
"CHECKBOX"			{ COUNT; return(CHECKBOX); }
"COMBOBOX"			{ COUNT; return(COMBOBOX); }
"CONTROL"			{ COUNT; return(CONTROL); }
"CTEXT"				{ COUNT; return(CTEXT); }
"DEFPUSHBUTTON"		{ COUNT; return(DEFPUSHBUTTON); }
"EDITTEXT"			{ COUNT; return(EDITTEXT); }
"GROUPBOX"			{ COUNT; return(GROUPBOX); }
"LISTBOX"			{ COUNT; return(LISTBOX); }
"LTEXT"				{ COUNT; return(LTEXT); }
"PUSHBOX"			{ COUNT; return(PUSHBOX); }
"PUSHBUTTON"		{ COUNT; return(PUSHBUTTON); }
"RADIOBUTTON"		{ COUNT; return(RADIOBUTTON); }
"RTEXT"				{ COUNT; return(RTEXT); }
"SCROLLBAR"			{ COUNT; return(SCROLLBAR); }
"STATE3"			{ COUNT; return(STATE3); }

"POPUP"				{ COUNT; return(POPUP); }
"MENUITEM"			{ COUNT; return(MENUITEM); }
"BUTTON"			{ COUNT; return(BUTTON); }
"SEPARATOR"			{ COUNT; return(SEPARATOR); }
"MFT_SEPARATOR"		{ COUNT; return(MFTSEPARATOR); }
"DLGINIT"			{ COUNT; return(DLGINIT); }

"BLOCK"				{ COUNT; return(BLOCK); }
"VALUE"				{ COUNT; return(VALUE); }

"MOVEABLE"			{ COUNT; return(MOVEABLE_FLAG); }
"FIXED"				{ COUNT; return(FIXED_FLAG); }
"PURE"				{ COUNT; return(PURE_FLAG); }
"IMPURE"			{ COUNT; return(IMPURE_FLAG); }
"PRELOAD"			{ COUNT; return(PRELOAD_FLAG); }
"LOADONCALL"		{ COUNT; return(LOADONCALL_FLAG); }
"DISCARDABLE"		{ COUNT; return(DISCARDABLE_FLAG); }

"NOT"				{ COUNT; return(NOT_OP); }

	/*
	  =========================
	         identifiers
	  =========================
	*/
{L}({L}|{D})*	{ COUNT; return(IDENTIFIER); }

	/*
	  =========================
	           numbers
	  =========================
	*/
[+-]?0[xX]{H}+{IS}?	{ COUNT; return(CONSTANT); }
[+-]?0{D}+{IS}?		{ COUNT; return(CONSTANT); }
[+-]?{D}+{IS}?		{ COUNT; return(CONSTANT); }

	/*
	  =========================
	           strings
	  =========================
	*/
\"(\\.|\"\"|[^\\"])*\"	{ COUNT; return(STRING); }

	/*
	  =========================
	          operators
	  =========================
	*/
">>"			{ COUNT; return(RIGHT_OP); }
"<<"			{ COUNT; return(LEFT_OP); }
"&&"			{ COUNT; return(AND_OP); }
"||"			{ COUNT; return(OR_OP); }
"<="			{ COUNT; return(LE_OP); }
">="			{ COUNT; return(GE_OP); }
"=="			{ COUNT; return(EQ_OP); }
"!="			{ COUNT; return(NE_OP); }
","				{ COUNT; return(','); }
"."				{ COUNT; return('.'); }
"="				{ COUNT; return('='); }
"("				{ COUNT; return('('); }
")"				{ COUNT; return(')'); }
"["				{ COUNT; return('['); }
"]"				{ COUNT; return(']'); }
"&"				{ COUNT; return('&'); }
"!"				{ COUNT; return('!'); }
"~"				{ COUNT; return('~'); }
"-"				{ COUNT; return('-'); }
"+"				{ COUNT; return('+'); }
"*"				{ COUNT; return('*'); }
"/"				{ COUNT; return('/'); }
"%"				{ COUNT; return('%'); }
"<"				{ COUNT; return('<'); }
">"				{ COUNT; return('>'); }
"^"				{ COUNT; return('^'); }
"|"				{ COUNT; return('|'); }
"$"				{ COUNT; return('$'); }

[ \t\v\n\f]		{ COUNT; }
.				{ return SCANERR; }

%%

int yywrap()
{
	return(1);
}

void count(yy_parse_state *state)
{
	int i;

	for (i = 0; yytext[i] != '\0'; i++) {
		if (yytext[i] == '\n') {
			state->curColumn = 1;
			state->curLine++;
		} else if (yytext[i] == '\t') {
			state->curColumn += 8 - ((state->curColumn-1) % 8);
		} else {
			state->curColumn++;
		}
	}
}

int yylex(int *yylval, void *yyinput)
{
	yy_parse_state *state = (yy_parse_state*)yyinput;
	int directive;
	int token;

	for (;;) {
		token = _yylex(yylval, yyinput);
		if (token == ELIF_DIRECTIVE
		 || token == IF_DIRECTIVE
		 || token == IFDEF_DIRECTIVE
		 || token == IFNDEF_DIRECTIVE) 
		{
			state->preprocessorParams.clear();
			directive = token;
			continue;
		} 
		if (token == DEFINE_DIRECTIVE
		 || token == ELSE_DIRECTIVE
		 || token == ENDIF_DIRECTIVE
		 || token == INCLUDE_DIRECTIVE
		 || token == PRAGMA_DIRECTIVE
		 || token == UNDEF_DIRECTIVE) 
		{
			state->preprocessorParams.clear();
			directive = token;
			token = DIRECTIVE_END;
		} 
		if (token == DIRECTIVE_PARAM) {
			state->preprocessorParams += yytext;
		} else if (token == DIRECTIVE_END) {
			if (directive == ELIF_DIRECTIVE) {
				if (state->conditionStack.size() > 1) {
					int cond = state->conditionStack.top(); 
					state->conditionStack.pop();
					if (cond == 0 && state->conditionStack.top() != 0)
						state->trigSkip = -1;
					cond = checkCondition(state->preprocessorParams.c_str(), true);
					if (state->conditionStack.top() != 0 && cond == 0)
						state->trigSkip = 1;
					state->conditionStack.push(cond);
					state->conditionState = cond;
				} else {
					yyerror("unmatched #elif");
					return -1;
				}
			} else if (directive == ELSE_DIRECTIVE) {
				if (state->conditionStack.size() > 1) {
					int cond = state->conditionStack.top(); 
					state->conditionStack.pop();
					if (cond == 0 && state->conditionStack.top() != 0)
						state->trigSkip = -1;
					cond = ((cond == 1) ? 0 : ((cond == 0) ? 1 : -1));
					if (state->conditionStack.top() != 0 && cond == 0)
						state->trigSkip = 1;
					state->conditionStack.push(cond);
					state->conditionState = cond;
				} else {
					yyerror("unmatched #else");
					return -1;
				}
			} else if (directive == ENDIF_DIRECTIVE) {
				if (state->conditionStack.size() > 1) {
					int cond = state->conditionStack.top(); 
					state->conditionStack.pop();
					if (cond == 0 && state->conditionStack.top() != 0)
						state->trigSkip = -1;
				} else {
					yyerror("unmatched #endif");
					return -1;
				}
			} else if (directive == IF_DIRECTIVE) {
				int cond = checkCondition(state->preprocessorParams.c_str(), true);
				if (state->conditionStack.top() != 0 && cond == 0)
					state->trigSkip = 1;
				state->conditionStack.push(cond);
				state->conditionState = cond;
			} else if (directive == IFDEF_DIRECTIVE) {
				int cond = checkCondition(state->preprocessorParams.c_str(), false);
				if (state->conditionStack.top() != 0 && cond == 0)
					state->trigSkip = 1;
				state->conditionStack.push(cond);
				state->conditionState = cond;
			} else if (directive == IFNDEF_DIRECTIVE) {
				int cond = checkCondition(state->preprocessorParams.c_str(), false);
				cond = ((cond == 1) ? 0 : ((cond == 0) ? 1 : -1));
				if (state->conditionStack.top() != 0 && cond == 0)
					state->trigSkip = 1;
				state->conditionStack.push(cond);
				state->conditionState = cond;
			}
		} else {
			break;
		}
	}

	return token;
}
