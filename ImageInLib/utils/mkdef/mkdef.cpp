#include <Windows.h>
#include <io.h>
#include <wchar.h>

#include <stdio.h>
#include <string.h>
#include <direct.h>
#include <stdlib.h>
#include <io.h>
#include <stdio.h>
#include <conio.h>
#include <direct.h>
#include <stdlib.h>
#include <ctype.h>
#include <dos.h>
#include <malloc.h>
#include <ImageHlp.h>

typedef void* (__cdecl *malloc_func_t)(size_t);
typedef void  (__cdecl *free_func_t)(void*);
extern "C" char* __cdecl __unDName(char* buffer, const char* mangled, int buflen, malloc_func_t memget, free_func_t memfree, unsigned short int flags);

#define LINE_LENGTH 1024

int		verbose = FALSE;
int 	export_cdecl = TRUE;
int 	ordinals_only = FALSE;

char	*defFile = NULL;
char	*mapFile = NULL;

typedef enum _symtype
{
	A_UNKNOWN,
	A_C_DECL,
	A_DATA,
	A_FUNC,
} SYMTYPE;

typedef struct _symbol
{
	char			*name;
	char			*undname;
	SYMTYPE			type;
	struct _symbol	*next;
} SYMBOL;

void getSymbolsFromMapFile(char	*mapFile,
						   int		verbose);

int getNumSymbols();

unsigned long addSymbol(char *name, SYMTYPE type, char *undname);

SYMBOL *getNextSymbol(SYMBOL	*sym);

void removeSymbols();

static unsigned long	numSymbols = 0;
static SYMBOL			*firstSymbol = NULL;
static SYMBOL			*lastSymbol = NULL;

int getNumSymbols()
{
	return numSymbols;
}

unsigned long addSymbol(char *name,	SYMTYPE	type, char *undname)
{
	SYMBOL	*s;
	
	s = (SYMBOL *)malloc(sizeof(SYMBOL));
	if (s == NULL)
		return FALSE;
		
	s->name = _strdup(name);
	s->undname = undname ? _strdup(undname) : NULL;
	s->type = type;
	s->next = NULL;
	
	if (lastSymbol != NULL)
		lastSymbol->next = s;
	lastSymbol = s;
	if (firstSymbol == NULL)
		firstSymbol = s;
	numSymbols++;

	return numSymbols;
}

SYMBOL *getNextSymbol(SYMBOL *sym)
{
	if (sym == NULL)
		return firstSymbol;
	return sym->next;
}

void removeSymbols()
{
	if (firstSymbol != NULL)
	{
		SYMBOL	*sym = NULL;
		SYMBOL	*symD = NULL;
		
		sym = firstSymbol;
		while(TRUE)
		{
			symD = sym;
			sym = getNextSymbol(sym);
			if (symD != NULL)
			{
				if (symD->name != NULL)
					free(symD->name);
				if (symD->undname != NULL)
					free(symD->undname);
				free(symD);
			}
			else
				break;
		}
		
		firstSymbol = NULL;
		lastSymbol = NULL;
	}
}

void getSymbolsFromMapFile(char	*mapFile,
						   int		verbose)
{
	FILE	*fp;
	char	line[LINE_LENGTH];
	char	*p;
	
	// open the map file
	
	fp = fopen(mapFile, "r");
	if (fp != NULL)
	{
		int	matchLine = FALSE;

		// search for " Address   Export                  Alias"
		// search for "  Address         Publics by Value              Rva+Base       Lib:Object"
	
		if (verbose)
			printf("Searching for start of exported names\n");
		while(TRUE)
		{
			p = fgets(line, (LINE_LENGTH - 1), fp);
			if (p == NULL)
				break;
		
			if (strncmp(line, "  Address", 9) == 0)
			{
				// check remainder of line has:-
				//	Publics by Value              
				//	Rva+Base     
				//	Lib:Object
				// all separated by spaces

				char	*tp;

				tp = strstr(line, "Publics by Value");
				if (tp != NULL)
				{
					tp = strstr(tp, "Rva+Base");
					if (tp != NULL)
					{
						tp = strstr(tp, "Lib:Object");
						if (tp != NULL)
						{
							matchLine = TRUE;
						}
					}
				}
			}

			if (matchLine)
			{
				// skip blank line
				
				p = fgets(line, (LINE_LENGTH - 1), fp);
				break;
			}
		}
		
		if (p != NULL)
		{
			char	line1[LINE_LENGTH];
			
			// load lines, 1 line at a time.
			// each line has the following format:-
			//
			// 0001:00000000       _DllMain@12                10001000 f UTIL.OBJ
			//
			// space, address, spaces, decorated-name, spaces, address, space, flagFunction, space, flagInline, space, obj or module (dll) name
			// flagFunction specifies f for function (I think - data doesn't have it)
			// flagInline specifies i for inline
			// the decorated name is exported unless it has the prefix __imp__
			
			if (verbose)
				printf("Parsing export names\n");
			while(TRUE)
			{
				// fetch lines if required
				// or shuffle them
				
				//AfxCheckMemory();
				p = fgets(line1, (LINE_LENGTH - 1), fp);
				if (p == NULL)
					break;

				// is this the correct line?
				
				// data we require is on first line
				// check that line starts with a space and a pointer address expressed as two hex words with a colon
				// " xxxx:xxxx",followed by one space, a name, an alias (often the same as the name)
				// the name and the alias may be a C name or a C++ name, we are only interested in the C names here
					
				if (line1[0] == ' ' && 
					isxdigit(line1[1]) &&
					isxdigit(line1[2]) &&
					isxdigit(line1[3]) &&
					isxdigit(line1[4]) &&
					line1[5] == ':' &&
					isxdigit(line1[6]) &&
					isxdigit(line1[7]) &&
					isxdigit(line1[8]) &&
					isxdigit(line1[9]) &&
					isxdigit(line1[10]) &&
					isxdigit(line1[11]) &&
					isxdigit(line1[12]) &&
					isxdigit(line1[13]) &&
					isspace(line1[14]))
				{
					// first name is separated from second name by whitespace
					
					char	*name = NULL;
					char    *undname = NULL;
					char	*address = NULL;
					char	*flagFunction = NULL;
					char	*flagInline = NULL;
					char	*module = NULL;
					char	*pq;
					int		isWep, isCPP, isOBJ, isDLL, isDestructor, isRTTI, isVFTable, isString, isImport, isFunction;
						
					// skip spaces

					name = &line1[14];
					while(isspace(*name))
						name++;
					if (*name == '\0')
						break;
						
					// terminate the name

					pq = strchr(name, ' ');
					if (*pq == '\0')
						break;
					*pq = '\0';
					
					// get the address of the name

					address = pq + 1;
					while(isspace(*address))
						address++;
					if (*address == '\0')
						break;

					// terminate the address

					pq = strchr(address, ' ');
					if (*pq == '\0')
						break;
					*pq = '\0';

					// get the flagFunction

					flagFunction = pq + 1;
					*(flagFunction + 1) = '\0';

					// get the flagInline

					flagInline = pq + 3;
					*(flagInline + 1) = '\0';

					// get the object/module name

					module = flagInline + 2;

					isDLL = ((strstr(module, ".DLL") != NULL) || (strstr(module, ".dll") != NULL));
					isOBJ = ((strstr(module, ".OBJ") != NULL) || (strstr(module, ".obj") != NULL));
					if (isOBJ)
					{
						// check that this object is not really in a DLL

						if (strstr(module, ":") != NULL)
						{
							isOBJ = FALSE;
							isDLL = TRUE;
						}
					}

					isDestructor = ((strncmp(name, "??_E", 4) == 0) ||	// is ??_E destructor
									(strncmp(name, "??_G", 4) == 0));	// is ??_G scalar destructor
					isRTTI = (strncmp(name, "??_R", 4) == 0);
					isVFTable = (strncmp(name, "??_7", 4) == 0);
					isString = (strncmp(name, "??_C", 4) == 0);
					isImport = (strncmp(name, "__imp__", 7) == 0);

					// check that name is not 'WEP' or a C++ name
					isWep = (strcmp(name, "WEP") == 0);					
					isCPP = name[0] == '?';
					isFunction = flagFunction && *flagFunction == 'f';

					if (verbose)
						printf("Potential C export: %s\n", name);

					if (!isWep && isOBJ && !isImport)
					{
						if (!isCPP)
						{
							if ((strncmp(name, "_ExtRawDllMain", 14) != 0) &&
								(strncmp(name, "_DllMain", 8) != 0))
							{
								if (verbose)
								{
									printf("C export: %s\n", name);
									printf("  address: %s\n", address);
									printf("  flag(Function == f): %s\n", flagFunction);
									printf("  flag(Inline == i): %s\n", flagInline);
									printf("  module: %s\n\n", module);
								}

								addSymbol(name, A_C_DECL, NULL);
							}
						}
						else if (isCPP)
						{
							if (!isDestructor && !isString)
							{
								if (verbose)
								{
									printf("C++ export: %s\n", name);
									printf("  address: %s\n", address);
									printf("  flag(Function == f): %s\n", flagFunction);
									printf("  flag(Inline == i): %s\n", flagInline);
									printf("  module: %s\n\n", module);
								}

								undname = __unDName(NULL, name, 0, malloc, free, 0);
								addSymbol(name, isFunction ? A_FUNC : A_DATA, undname);
							}
						}
					}
				}
				else
					break;	// drop out of export parsing
			}
		}
		else
		{
			if (verbose)
				printf("No C or C++ names found!\n");
		}

		// close map file
	
		if (verbose)
			printf("Parsing C and C++ decorated names complete\n");
		fclose(fp);
	}
	else
	{
		if (verbose)
			printf("Failed to open map file %s\n", mapFile);
	}
}

void outputToDefFile(char	*defFile)
{
	FILE			*fp;
	FILE			*fp2;
	unsigned long	id = 100;
	SYMBOL			*sym = NULL;
	char			*p;
	char			line[LINE_LENGTH];
	char			tempFile[] = "fnXXXXXX";
		
	// open the def file
	
	//tempFile = 
		_mktemp(&tempFile[0]);
	
	fp = fopen(defFile, "r");
	fp2 = fopen(tempFile, "w+");
	if (fp2 != NULL)
	{
		if (fp != NULL)
		{
			if (verbose)
				printf("Parsing definition file\n");
			while(TRUE)
			{
				p = fgets(line, (LINE_LENGTH - 1), fp);
				if (p == NULL)
					break;
            
				if (strncmp(line, "EXPORTS", 7) == 0)
				{
					fputs(line, fp2);
			        	
					break;
				}
						
				fputs(line, fp2);
			}
		}
		else
		{
			if (fp == NULL && verbose)
				printf("Failed to open definitions file %s\n", defFile);
			p = NULL;
		}
		
		if (p == NULL)
		{
			fputs("EXPORTS\n", fp2);
			fputs("; mkdef-generated-exports-start\n", fp2);
		}
		else
		{
			while(TRUE)
			{
				p = fgets(line, (LINE_LENGTH - 1), fp);
				if (p == NULL)
					break;

				if (strncmp(line, "; mkdef-generated-exports-start", strlen("; mkdef-generated-exports-start")) == 0)
				{
					fputs(line, fp2);
					break;
				}
				fputs(line, fp2);
			}			

			if (p == NULL)
				fputs("; mkdef-generated-exports-start", fp2);
		}	

		// write all symbols that have been found in map, and have a class name in the class list

		if (verbose)
			printf("Writing definition file\n");
		
		sym = NULL;
		while(TRUE)
		{
			int	ok = FALSE;
			
			sym = getNextSymbol(sym);
			if (sym == NULL)
				break;

			// not a symbol in a C++ class
			// check if a C declaration
				
			if (sym->type == A_C_DECL)
				ok = export_cdecl;
			else if (sym->type == A_DATA)
				ok = TRUE;
			else if (sym->type == A_FUNC)
				ok = TRUE;
			
			if (ok)
			{
				char	t[100];

				sprintf(t, " @ %ld", id);
				fputs("    ", fp2);
				fputs(sym->name, fp2);
				fputs(t, fp2);
				if (ordinals_only)
					fputs(" NONAME", fp2);
				if (sym->type == A_DATA)
					fputs(" DATA", fp2);
				if (sym->undname)
				{
					fputs("     ; ", fp2);
					fputs(sym->undname, fp2);
				}
				fputs("\n", fp2);
				id++;
					
				if (verbose)
					printf("        %s%s", sym->name, t);
			}
		}

		// add terminator comment
		fputs("; mkdef-generated-exports-end\n", fp2);

		// continue copying data
		
		if (fp != NULL)
		{
			// find end of EXPORT section added by mkdef
			while(TRUE)
			{
				p = fgets(line, (LINE_LENGTH - 1), fp);
				if (p == NULL)
					break;

				if (strncmp(line, "; mkdef-generated-exports-end", strlen("; mkdef-generated-exports-end")) == 0)
					break;

                // or look for another directive...
                
				if (strncmp(line, "EXPORTS", 7) == 0)
					break;
			}

			// and copy the rest of the file
			
			if (p != NULL)
			{
				while(TRUE)
				{
					p = fgets(line, (LINE_LENGTH - 1), fp);
					if (p == NULL)
						break;
                    
                    fputs(line, fp2);
				}
			}
			
			fclose(fp);
        }
        
		// close def file
	
		if (verbose)
			printf("Definition file output complete\n");
		fclose(fp2);
		
		// now copy tempFile onto defFile
		
		remove(defFile);
		rename(tempFile, defFile);
	}
	
	if (fp2 == NULL && verbose)
		printf("Failed to open temporary definitions file %s\n", tempFile);
}


int main(int	argc,
		 char	*argv[])
{
	int		needHelp = FALSE;
	int		i;
	int		clean = FALSE;
	
	// handle options
	
	verbose = FALSE;
	if (argc > 1)
	{
		for(i = 1; i < argc;)
		{
			int	concat = 1;

			if (_stricmp(argv[i], "--verbose") == 0)
				verbose = TRUE;
			else if (_stricmp(argv[i], "-v") == 0)
				verbose = TRUE;
			else if (_stricmp(argv[i], "--help") == 0)
				needHelp = TRUE;
			else if (_stricmp(argv[i], "--cdecl") == 0)
				export_cdecl = TRUE;
			else if (_stricmp(argv[i], "--no-cdecl") == 0)
				export_cdecl = FALSE;
			else if (_stricmp(argv[i], "--ordinals-only") == 0)
				ordinals_only = TRUE;
			else if (_stricmp(argv[i], "--clean") == 0)
				clean = TRUE;
			else
				concat = 0;
				
			if (concat == 0)
				i++;
			else
			{
				while(concat > 0)
				{
					for(int k = i; k < argc - 1; k++)
						argv[k] = argv[k + 1];
					argc--;
					concat--;
				}
			}
		}
	}
	
	// handle help message
	
	if (!clean && (argc <= 2 ||	needHelp))
	{
		printf(".MAP file to .DEF file converter\n");
		printf("Usage: %s [options] [file.def] [file.map]\n", argv[0]);
		printf("Options are:-\n" \
			   "--verbose|-v\tdisplay verbose output\n" \
			   "--clean\t\tclean definition file and quit\n" \
			   "--help\t\tdisplay help message\n" \
			   "--ordinals-only\texport by ordinals only\n" \
			   "--cdecl\texport C functions\n" \
			   "--no-cdecl\tdo not export C functions\n");
		return 0;
	}

	defFile = _strdup(argv[1]);
	mapFile = _strdup(argv[2]);

	if (verbose)
	{		
		printf("Definitions file:%s\n", defFile);
		printf("Map file:%s\n\n", mapFile);
    }
    
	// do work

	if (!clean)
	{	
		getSymbolsFromMapFile(mapFile, verbose);

		if (verbose)
			printf("Total number of symbols:%d\n\n", getNumSymbols());
	}
	outputToDefFile(defFile);
	
	// tidyup
	
	removeSymbols();
	
	if (mapFile != NULL)
		free(mapFile);
	if (defFile != NULL)
		free(defFile);
		
	return 0;
}
