/*
 *  X86 code generator for TCC
 * 
 *  Copyright (c) 2001-2004 Fabrice Bellard
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include "tcc.h"
#include <assert.h>
#include "svm-opcodes.h"


#ifdef TCC_TARGET_SVM32
const int PTR_SIZE = 4;
const int EM_TCC_TARGET  = EM_386;
const int R_PC_32        = R_386_PC32;
const int R_DATA_32      = R_386_32;
const int R_DATA_PTR     = R_386_32;
const int R_JMP_SLOT     = R_386_JMP_SLOT;
const int R_COPY         = R_386_COPY;
const int ELF_START_ADDR = 0x08048000;
const int ELF_PAGE_SIZE  = 0x1000;
#endif

#ifdef TCC_TARGET_SVM64
const int PTR_SIZE = 8;
const int EM_TCC_TARGET  = EM_X86_64;
const int R_PC_32        = R_X86_64_PC32;
const int R_DATA_32      = R_X86_64_32;
const int R_DATA_PTR     = R_X86_64_64;
const int R_JMP_SLOT     = R_X86_64_JUMP_SLOT;
const int R_COPY         = R_X86_64_COPY;
const int ELF_START_ADDR = 0x400000;
const int ELF_PAGE_SIZE  = 0x200000;
#endif


ST_DATA const int reg_classes[NB_REGS] = {
    /* R00 */ RC_INT | RC_FLOAT | RC_R00,
    /* R08 */ RC_INT | RC_FLOAT | RC_R08,
    /* R10 */ RC_INT | RC_FLOAT,
    /* R18 */ RC_INT | RC_FLOAT,
    /* R20 */ RC_INT | RC_FLOAT,
    /* R28 */ RC_INT | RC_FLOAT,
    /* R30 */ RC_INT | RC_FLOAT,
    /* R38 */ RC_INT | RC_FLOAT,
    /* R40 */ RC_INT | RC_FLOAT,
    /* R48 */ RC_INT | RC_FLOAT,
    /* R50 */ RC_INT | RC_FLOAT,
    /* R58 */ RC_INT | RC_FLOAT,
    /* R60 */ RC_INT | RC_FLOAT,
    /* R68 */ RC_INT | RC_FLOAT,
    /* R70 */ RC_INT | RC_FLOAT,
    /* R78 */ RC_INT | RC_FLOAT,
    /* R80 */ RC_INT | RC_FLOAT,
    /* R88 */ RC_INT | RC_FLOAT,
    /* R90 */ RC_INT | RC_FLOAT,
    /* R98 */ RC_INT | RC_FLOAT,
    /* RA0 */ RC_INT | RC_FLOAT,
    /* RA8 */ RC_INT | RC_FLOAT,
    /* RB0 */ RC_INT | RC_FLOAT,
    /* RB8 */ RC_INT | RC_FLOAT,
    /* RC0 */ RC_INT | RC_FLOAT,
    /* RC8 */ RC_INT | RC_FLOAT,
    /* RD0 */ RC_INT | RC_FLOAT,
    /* RD8 */ RC_INT | RC_FLOAT,
    /* RE0 */ RC_INT | RC_FLOAT,
    /* RE8 */ RC_INT | RC_FLOAT,
    /* RF0 */ RC_INT | RC_FLOAT,
    /* RF8 */ RC_INT | RC_FLOAT,
};


static unsigned long func_sub_sp_offset;
static int func_ret_sub;

static int is64_type(int t)
{
    return ((t & VT_BTYPE) == VT_LLONG) ||
		   ((t & VT_BTYPE) == VT_DOUBLE) ||
		   ((t & VT_BTYPE) == VT_LDOUBLE) ||
		   ((PTR_SIZE == 8) && (((t & VT_BTYPE) == VT_PTR) || ((t & VT_BTYPE) == VT_FUNC)));
}

ST_FUNC void gen1(int c)
{
    int ind1;
    ind1 = ind + 1;
    if (ind1 > cur_text_section->data_allocated)
        section_realloc(cur_text_section, ind1);
    cur_text_section->data[ind] = c;
    ind = ind1;
}

ST_FUNC void gen2(int v)
{
    gen1(v);
    gen1(v >> 8);
}

ST_FUNC void gen3(int c)
{
    gen1(c);
    gen1(c >> 8);
    gen1(c >> 16);
}

ST_FUNC void gen4(int c)
{
    gen1(c);
    gen1(c >> 8);
    gen1(c >> 16);
    gen1(c >> 24);
}

ST_FUNC void gen8(int64_t c)
{
    gen1(c);
    gen1(c >> 8);
    gen1(c >> 16);
    gen1(c >> 24);
    gen1(c >> 32);
    gen1(c >> 40);
    gen1(c >> 48);
    gen1(c >> 56);
}

/* output constant with relocation if 'r & VT_SYM' is true */
ST_FUNC void gen_addr32(int r, Sym *sym, int c)
{
    if (r & VT_SYM)
        greloc(cur_text_section, sym, ind, R_DATA_32);
    gen4(c);
}

/* output constant with relocation if 'r & VT_SYM' is true */
ST_FUNC void gen_addr64(int r, Sym *sym, int64_t c)
{
    if (r & VT_SYM)
        greloc(cur_text_section, sym, ind, R_X86_64_64);
    gen8(c);
}

/* output constant with relocation if 'r & VT_SYM' is true */
ST_FUNC void gen_addrpc32(int r, Sym *sym, int c)
{
    if (r & VT_SYM)
        greloc(cur_text_section, sym, ind, R_PC_32);
    gen4(c - 4);
}

/* output a symbol and patch all calls to it */
ST_FUNC void gsym_addr(int t, int a)
{
    int n, *ptr;
    while (t)
	{
        ptr = (int *)(cur_text_section->data + t);
		if (ind == a)
			tcc_listing_print("@@%x:\n", t, cur_text_section->name);
		else
			tcc_listing_print("@@%x=%s+0x%05x\n", t, cur_text_section->name, a);
        n = *ptr; /* next value */
        *ptr = a - t - 4;
        t = n;
    }
}

ST_FUNC void gsym(int t)
{
    gsym_addr(t, ind);
}

/* psym is used to put a data field which is a reference to a symbol.
   return the address of the data field */
ST_FUNC int psym(int s)
{
    int ind1;
    ind1 = ind + 4;
    if (ind1 > cur_text_section->data_allocated)
        section_realloc(cur_text_section, ind1);
    *(int *)(cur_text_section->data + ind) = s;
    s = ind;
    ind = ind1;
    return s;
}

static void gen_load_const(int r, CValue *val, int type)
{
	int ll = is64_type(type);
	r = REG_VALUE(r)*8;
	if (!ll)
	{
		if (val->ul == (char)val->ul)
		{
			gen1(SVMOPCODE_MOVL_IM8S_R);
			gen1(val->ul);
			gen1(r);
			tcc_listing_code(3, -1, "movl $%d, %%r%02X", (int)val->ul, r);
		}
		else if (val->ul == (short)val->ul)
		{
			gen1(SVMOPCODE_MOVL_IM16S_R);
			gen2(val->ul);
			gen1(r);
			tcc_listing_code(4, -1, "movl $%d, %%r%02X", (int)val->ul, r);
		}
		else
		{
			gen1(SVMOPCODE_MOVL_IMM_R);
			gen4(val->ul);
			gen1(r);
			tcc_listing_code(6, -1, "movl $%d, %%r%02X", (int)val->ul, r);
		}
	}
	else
	{
		gen1(SVMOPCODE_MOVQ_IMM_R);
		gen8(val->ull);
		gen1(r);
		tcc_listing_code(10, -1, "movq $%lld, %%r%02X", val->ull, r);
	}
}

ST_FUNC void gen_load_reg(int dst, int src, int type)
{
	int op;
	char *mne;
	int ll = is64_type(type);
	if (is64_type(type))
	{
		op = SVMOPCODE_MOVQ_R_R;
		mne = "movq";
	}
	else if (dst & TREG_CAST64)
	{
		 if ((type & (VT_BTYPE | VT_UNSIGNED)) == VT_INT)
		 {
			op = SVMOPCODE_MOVSLQ_R_R;
			mne = "movslq";
		 }
		 else
		 {
			op = SVMOPCODE_MOVZLQ_R_R;
			mne = "movzlq";
		 }
	}
	else
	{
		op = SVMOPCODE_MOVL_R_R;
		mne = "movl";
	}
	dst = REG_VALUE(dst)*8;
	src = REG_VALUE(src)*8;
	gen1(op);
	gen1(src);
	gen1(dst);
	tcc_listing_code(3, -1, "%s %%r%02X, %%r%02X",mne, src, dst);
}

static void gen_load_ind(int dst, int src, int type)
{
	int op;
	char *mne;
	dst = REG_VALUE(dst)*8;
	src = REG_VALUE(src)*8;
	if (is64_type(type))
	{
		op = SVMOPCODE_MOVQ_IND_R;
		mne = "movq";
	}
	else if ((type & VT_TYPE) == VT_BYTE || (type & VT_TYPE) == VT_BOOL)
	{
		op = SVMOPCODE_MOVSBL_IND_R;
		mne = "movsbl";
	}
	else if ((type & VT_TYPE) == (VT_BYTE | VT_UNSIGNED))
	{
		op = SVMOPCODE_MOVZBL_IND_R;
		mne = "movzbl";
	}
	else if ((type & VT_TYPE) == VT_SHORT)
	{
		op = SVMOPCODE_MOVSWL_IND_R;
		mne = "movswl";
	}
	else if ((type & VT_TYPE) == (VT_SHORT | VT_UNSIGNED))
	{
		op = SVMOPCODE_MOVZWL_IND_R;
		mne = "movzwl";
	}
	else
	{
		op = SVMOPCODE_MOVL_IND_R;
		mne = "movl";
	}
	gen1(op);
	gen1(src);
	gen1(dst);
	tcc_listing_code(3, -1, "%s (%%r%02X), %%r%02X", mne, src, dst);
}

static void gen_load_rip(int r, int value, Sym *sym, int disp, int type)
{
	int op;
	char *mne;
	r = REG_VALUE(r)*8;
	if (is64_type(type))
	{
		op = SVMOPCODE_MOVQ_RIP_R;
		mne = "movq";
	}
	else if ((type & VT_TYPE) == VT_BYTE || (type & VT_TYPE) == VT_BOOL)
	{
		op = SVMOPCODE_MOVSBL_RIP_R;
		mne = "movsbl";
	}
	else if ((type & VT_TYPE) == (VT_BYTE | VT_UNSIGNED))
	{
		op = SVMOPCODE_MOVZBL_RIP_R;
		mne = "movzbl";
	}
	else if ((type & VT_TYPE) == VT_SHORT)
	{
		op = SVMOPCODE_MOVSWL_RIP_R;
		mne = "movswl";
	}
	else if ((type & VT_TYPE) == (VT_SHORT | VT_UNSIGNED))
	{
		op = SVMOPCODE_MOVZWL_RIP_R;
		mne = "movzwl";
	}
	else
	{
		op = SVMOPCODE_MOVL_RIP_R;
		mne = "movl";
	}
	gen1(op);
    gen_addrpc32(value, sym, disp);
	gen1(r);
	if (value & VT_SYM)
	{
	    char *symname = get_tok_str(sym->v, NULL);
		if (disp)
			tcc_listing_code(6, 1, "%s %s+0x%x(%%rip), %%r%02X", mne, symname, disp, r);
		else
			tcc_listing_code(6, 1, "%s %s(%%rip), %%r%02X", mne, symname, r);
	}
	else
	{
		tcc_listing_code(6, -1, "%s 0x%x(%%rip), %%r%02X", mne, disp, r);
	}
}

static void gen_load_rbp(int r, int disp, int type)
{
	int op;
	char *mne;
	int disp8 = (disp == (char)disp);
	r = REG_VALUE(r)*8;
	if (is64_type(type))
	{
		op = disp8 ? SVMOPCODE_MOVQ_RBP8_R : SVMOPCODE_MOVQ_RBP_R;
		mne = "movq";
	}
	else if ((type & VT_TYPE) == VT_BYTE || (type & VT_BTYPE) == VT_BOOL)
	{
		disp8 = 0;
		op = SVMOPCODE_MOVSBL_RBP_R;
		mne = "movsbl";
	}
	else if ((type & VT_TYPE) == (VT_BYTE | VT_UNSIGNED))
	{
		disp8 = 0;
		op = SVMOPCODE_MOVZBL_RBP_R;
		mne = "movzbl";
	}
	else if ((type & VT_TYPE) == VT_SHORT)
	{
		disp8 = 0;
		op = SVMOPCODE_MOVSWL_RBP_R;
		mne = "movswl";
	}
	else if ((type & VT_TYPE) == (VT_SHORT | VT_UNSIGNED))
	{
		disp8 = 0;
		op = SVMOPCODE_MOVZWL_RBP_R;
		mne = "movzwl";
	}
	else
	{
		op = disp8 ? SVMOPCODE_MOVL_RBP8_R : SVMOPCODE_MOVL_RBP_R;
		mne = "movl";
	}
	gen1(op);
	if (disp8)
		gen1((char)disp);
	else
		gen4(disp);
	gen1(r);
	tcc_listing_code(disp8 ? 3 : 6, -1, "%s %d(%%rbp), %%r%02X", mne, disp, r);
}

static void gen_lea_rip(int dst, int src, Sym *sym, int disp)
{
	dst = REG_VALUE(dst)*8;
	gen1(SVMOPCODE_LEA_RIP_R);
    gen_addrpc32(src, sym, disp);
	gen1(dst);
	if (src & VT_SYM)
	{
	    char *symname = get_tok_str(sym->v, NULL);
		if (disp)
			tcc_listing_code(6, 1, "lea %s+0x%x(%%rip), %%r%02X", symname, disp, dst);
		else
			tcc_listing_code(6, 1, "lea %s(%%rip), %%r%02X", symname, dst);
	}
	else
	{
		tcc_listing_code(6, -1, "lea 0x%x(%%rip), %%r%02X", disp, dst);
	}
}

static void gen_lea_rbp(int dst, int src, Sym *sym, int disp)
{
	dst = REG_VALUE(dst)*8;
	gen1(SVMOPCODE_LEA_RBP_R);
    gen4(disp);
	gen1(dst);
	if (src & VT_SYM)
	{
	    char *symname = get_tok_str(sym->v, NULL);
		if (disp)
			tcc_listing_code(6, 1, "lea %s+0x%x(%%rbp), %%r%02X", symname, disp, dst);
		else
			tcc_listing_code(6, 1, "lea %s(%%rbp), %%r%02X", symname, dst);
	}
	else
	{
		tcc_listing_code(6, -1, "lea 0x%x(%%rbp), %%r%02X", disp, dst);
	}
}

static void gen_store_ind(int src, int dst, int type)
{
	int op;
	char *mne;
	dst = REG_VALUE(dst)*8;
	src = REG_VALUE(src)*8;
	if (is64_type(type))
	{
		op = SVMOPCODE_MOVQ_R_IND;
		mne = "movq";
	}
	else if ((type & VT_BTYPE) == VT_BYTE || (type & VT_BTYPE) == VT_BOOL)
	{
		op = SVMOPCODE_MOVB_R_IND;
		mne = "movb";
	}
	else if ((type & VT_BTYPE) == VT_SHORT)
	{
		op = SVMOPCODE_MOVW_R_IND;
		mne = "movw";
	}
	else
	{
		op = SVMOPCODE_MOVL_R_IND;
		mne = "movl";
	}
	gen1(op);
	gen1(src);
	gen1(dst);
	tcc_listing_code(3, -1, "%s %%r%02X, (%%r%02X)", mne, src, dst);
}

static void gen_store_rip(int r, int value, Sym *sym, int disp, int type)
{
	int op;
	char *mne;
	r = REG_VALUE(r)*8;
	if (is64_type(type))
	{
		op = SVMOPCODE_MOVQ_R_RIP;
		mne = "movq";
	}
	else if ((type & VT_BTYPE) == VT_BYTE || (type & VT_BTYPE) == VT_BOOL)
	{
		op = SVMOPCODE_MOVB_R_RIP;
		mne = "movb";
	}
	else if ((type & VT_BTYPE) == VT_SHORT)
	{
		op = SVMOPCODE_MOVW_R_RIP;
		mne = "movw";
	}
	else
	{
		op = SVMOPCODE_MOVL_R_RIP;
		mne = "movl";
	}
	gen1(op);
	gen1(r);
    gen_addrpc32(value, sym, disp);
	if (value & VT_SYM)
	{
	    char *symname = get_tok_str(sym->v, NULL);
		if (disp)
			tcc_listing_code(6, 1, "%s %%r%02X, %s+0x%x(%%rip)", mne, r, symname, disp);
		else
			tcc_listing_code(6, 1, "%s %%r%02X, %s(%%rip)", mne, r, symname);
	}
	else
	{
		tcc_listing_code(6, -1, "%s %%r%02X, 0x%x(%%rip)", mne, r, disp);
	}
}

static void gen_store_rbp(int r, int disp, int type)
{
	int op;
	char *mne;
	int disp8 = (disp == (char)disp);
	r = REG_VALUE(r)*8;
	if (is64_type(type))
	{
		op = disp8 ? SVMOPCODE_MOVQ_R_RBP8 : SVMOPCODE_MOVQ_R_RBP;
		mne = "movq";
	}
	else if ((type & VT_BTYPE) == VT_BYTE || (type & VT_BTYPE) == VT_BOOL)
	{
		op = disp8 ? SVMOPCODE_MOVB_R_RBP8 : SVMOPCODE_MOVB_R_RBP;
		mne = "movb";
	}
	else if ((type & VT_BTYPE) == VT_SHORT)
	{
		op = disp8 ? SVMOPCODE_MOVW_R_RBP8 : SVMOPCODE_MOVW_R_RBP;
		mne = "movw";
	}
	else
	{
		op = disp8 ? SVMOPCODE_MOVL_R_RBP8 : SVMOPCODE_MOVL_R_RBP;
		mne = "movl";
	}
	gen1(op);
	gen1(r);
	if (disp8)
		gen1((char)disp);
	else
		gen4(disp);
	tcc_listing_code(disp8 ? 3 : 6, -1, "%s %%r%02X, %d(%%rbp)", mne, r, disp);
}

static void gen_set_reg(int r, int cond)
{
	int op;
	char *mne;
	r = REG_VALUE(r)*8;
	switch (cond)
	{
	case TOK_ULT:
		op = SVMOPCODE_SETB;
		mne = "setb";
		break;
	case TOK_UGE:
		op = SVMOPCODE_SETAE;
		mne = "setae";
		break;
	case TOK_EQ:
		op = SVMOPCODE_SETE;
		mne = "sete";
		break;
	case TOK_NE:
		op = SVMOPCODE_SETNE;
		mne = "setne";
		break;
	case TOK_ULE:
		op = SVMOPCODE_SETBE;
		mne = "setbe";
		break;
	case TOK_UGT:
		op = SVMOPCODE_SETA;
		mne = "seta";
		break;
	case TOK_Nset:
		op = SVMOPCODE_SETS;
		mne = "sets";
		break;
	case TOK_Nclear:
		op = SVMOPCODE_SETNS;
		mne = "setns";
		break;
	case TOK_LT:
		op = SVMOPCODE_SETL;
		mne = "setl";
		break;
	case TOK_GE:
		op = SVMOPCODE_SETGE;
		mne = "setge";
		break;
	case TOK_LE:
		op = SVMOPCODE_SETLE;
		mne = "setle";
		break;
	case TOK_GT:
		op = SVMOPCODE_SETG;
		mne = "setg";
		break;
	case TOK_PE:
		op = SVMOPCODE_SETP;
		mne = "setp";
		break;
	case TOK_PO:
		op = SVMOPCODE_SETNP;
		mne = "setnp";
		break;
	}
	gen1(op);
	gen1(r);
	tcc_listing_code(2, -1, "%s %%r%02X", mne, r);
}

static int gen_branch_psym(int t, int cond)
{
	int ind1;
	int op;
	char *mne;
	switch (cond)
	{
	case TOK_ULT:
		op = SVMOPCODE_JB;
		mne = "jb";
		break;
	case TOK_UGE:
		op = SVMOPCODE_JAE;
		mne = "jae";
		break;
	case TOK_EQ:
		op = SVMOPCODE_JE;
		mne = "je";
		break;
	case TOK_NE:
		op = SVMOPCODE_JNE;
		mne = "jne";
		break;
	case TOK_ULE:
		op = SVMOPCODE_JBE;
		mne = "jbe";
		break;
	case TOK_UGT:
		op = SVMOPCODE_JA;
		mne = "ja";
		break;
	case TOK_Nset:
		op = SVMOPCODE_JS;
		mne = "js";
		break;
	case TOK_Nclear:
		op = SVMOPCODE_JNS;
		mne = "jns";
		break;
	case TOK_LT:
		op = SVMOPCODE_JL;
		mne = "jl";
		break;
	case TOK_GE:
		op = SVMOPCODE_JGE;
		mne = "jge";
		break;
	case TOK_LE:
		op = SVMOPCODE_JLE;
		mne = "jle";
		break;
	case TOK_GT:
		op = SVMOPCODE_JG;
		mne = "jg";
		break;
	case TOK_PE:
		op = SVMOPCODE_JP;
		mne = "jp";
		break;
	case TOK_PO:
		op = SVMOPCODE_JNP;
		mne = "jnp";
		break;
	}
	gen1(op);
	ind1 = ind;
	t = psym(t);
	tcc_listing_code(5, 1, "%s @@%x", mne, ind1);
	return t;
}

static void gen_jmp_addr(int a)
{
    int r;
    r = a - ind - 2;
    if (r == (char)r)
	{
		gen1(SVMOPCODE_JMPS);
		gen1(r);
		tcc_listing_code(2, -1, "jmp %s+0x%05x", cur_text_section->name, a);
	}
	else
	{
		gen1(SVMOPCODE_JMP);
		gen4(a - ind - 4);
		tcc_listing_code(5, -1, "jmp %s+0x%05x", cur_text_section->name, a);
	}
}

static int gen_jmp_psym(int t)
{
	int ind1;
	gen1(SVMOPCODE_JMP);
	ind1 = ind;
	t = psym(t);
	tcc_listing_code(5, 1, "jmp @@%x", ind1);
	return t;
}

static void gen_jmp_reg(int r)
{
	r = REG_VALUE(r)*8;
	gen1(SVMOPCODE_JMP_R);
	gen1(r);
	tcc_listing_code(2, -1, "jmp %%r%02X", r);
}

static void gen_call_reg(int r)
{
	r = REG_VALUE(r)*8;
	gen1(SVMOPCODE_CALL_R);
	gen1(r);
	tcc_listing_code(2, -1, "call %%r%02X", r);
}

static void gen_ncall_reg(int r, int args_size)
{
	r = REG_VALUE(r)*8;
	if (args_size == 0)
	{
		gen1(SVMOPCODE_NCALL0_R);
		gen1(r);
		tcc_listing_code(2, -1, "ncall0 %%r%02X", args_size, r);
	}
	else
	{
		gen1(SVMOPCODE_NCALL_IMM_R);
		gen2(args_size);
		gen1(r);
		tcc_listing_code(4, -1, "ncall $%d, %%r%02X", args_size, r);
	}
}

static void gen_ntcall_reg(int r, int args_size)
{
	r = REG_VALUE(r)*8;
	gen1(SVMOPCODE_NTCALL_IMM_R);
	gen2(args_size);
	gen1(r);
	tcc_listing_code(4, -1, "ntcall $%d, %%r%02X", args_size, r);
}

static void gen_jmp_rel(int r, Sym *sym, int64_t addr)
{
	gen1(SVMOPCODE_JMP);
	gen4(addr - 4);
	if (r & VT_SYM)
	{
	    char *symname = get_tok_str(sym->v, NULL);
        /* relocation case */
        greloc(cur_text_section, vtop->sym, 
			ind - 4, (PTR_SIZE==8) ? R_X86_64_PC32 : R_386_PC32);
		if (addr)
			tcc_listing_code(PTR_SIZE == 8 ? 10 : 6, -1, "jmp %s+0x%llx", symname, addr);
		else
			tcc_listing_code(PTR_SIZE == 8 ? 10 : 6, -1, "jmp %s", symname, addr);
	}
	else
	{
        /* put an empty PC32 relocation */
        put_elf_reloc(symtab_section, cur_text_section, 
			ind - 4, (PTR_SIZE==8) ? R_X86_64_PC32 : R_386_PC32, 0);
		tcc_listing_code(PTR_SIZE == 8 ? 10 : 6, -1, "jmp 0x%llx", addr);
	}
}

static void gen_call_rel(int r, Sym *sym, int64_t addr)
{
	gen1(SVMOPCODE_CALL);
	gen4(addr - 4);
	if (r & VT_SYM)
	{
	    char *symname = get_tok_str(sym->v, NULL);
        /* relocation case */
        greloc(cur_text_section, vtop->sym, 
			ind - 4, (PTR_SIZE == 8) ? R_X86_64_PC32 : R_386_PC32);
		if (addr)
			tcc_listing_code(5, 1, "call %s+0x%llx", symname, addr);
		else
			tcc_listing_code(5, 1, "call %s", symname, addr);
	}
	else
	{
        /* put an empty PC32 relocation */
        put_elf_reloc(symtab_section, cur_text_section, 
			ind - 4, (PTR_SIZE == 8) ? R_X86_64_PC32 : R_386_PC32, 0);
		tcc_listing_code(5, 1, "call 0x%llx", addr);
	}
}

static void gen_addsp(int val, int allowshort)
{
    if (allowshort && val == (char)val)
	{
		gen1(SVMOPCODE_ADDSP_IM8S);
        gen1(val);
		tcc_listing_code(2, -1, "add $%d, %%rsp", val);
    } else {
		gen1(SVMOPCODE_ADDSP_IMM);
        gen4(val);
		tcc_listing_code(5, -1, "add $%d, %%rsp", val);
    }
}

static void gen_subsp_reg(int r)
{
	r = REG_VALUE(r)*8;
	gen1(SVMOPCODE_SUBSP_R);
    gen1(r);
	tcc_listing_code(2, -1, "add %%r%02X, %%rsp", r);
}

static void gen_andsp_IM8S(int val)
{
	gen1(SVMOPCODE_ANDSP_IM8S);
    gen1(val);
	tcc_listing_code(2, -1, "and $%d, %%rsp", val);
}

static void gen_enter(int stacksize)
{
	gen1(SVMOPCODE_ENTER);
	if (stacksize == -1)
	{
		gen4(0);
		tcc_listing_code(5, 1, "enter %s@__LOCAL_SIZE__", funcname);
	}
	else
	{
		gen4(stacksize);
		tcc_listing_code(5, -1, "enter $%d", stacksize);
	}
}

static void gen_leave()
{
	gen1(SVMOPCODE_LEAVE);
	tcc_listing_code(1, -1, "leave");
}

static void gen_ret(int retsub)
{
	assert(retsub <= 65535);
	if (retsub)
	{
		gen1(SVMOPCODE_RETN);
		gen2(retsub);
		tcc_listing_code(3, -1, "retn 0x%x", retsub);
	}
	else
	{
		gen1(SVMOPCODE_RET);
		tcc_listing_code(1, -1, "ret");
	}
}

ST_FUNC void gen_ldsp_rbp(int disp)
{
	gen1(SVMOPCODE_LDSP_RBP);
	gen4(disp);
	tcc_listing_code(5, -1, "ldsp %d(%%rbp)", disp);
}

ST_FUNC void gen_stsp_rbp(int disp)
{
	gen1(SVMOPCODE_STSP_RBP);
	gen4(disp);
	tcc_listing_code(5, -1, "stsp %d(%%rbp)", disp);
}

ST_FUNC void gen_stsp_reg(int r)
{
	r = REG_VALUE(r)*8;
	gen1(SVMOPCODE_STSP_R);
	gen1(r);
	tcc_listing_code(2, -1, "stsp %%r%02X", r);
}

ST_FUNC void gen_push_reg(int r, int type)
{
	int op;
	char *mne;
	if (is64_type(type))
	{
		op = SVMOPCODE_PUSHQ_R;
		mne = "pushq";
	}
	else
	{
		op = SVMOPCODE_PUSHL_R;
		mne = "pushl";
	}
	r = REG_VALUE(r)*8;
	gen1(op);
	gen1(r);
	tcc_listing_code(2, -1, "%s %%r%02X", mne, r);
}

ST_FUNC void gen_push_reg32(int r, int type)
{
	assert(!is64_type(type));
	r = REG_VALUE(r)*8;
	gen1(SVMOPCODE_PUSH32_R);
	gen1(r);
	tcc_listing_code(2, -1, "push32 %%r%02X", r);
}

ST_FUNC void gen_imul_reg(int r, int multiplier, int type)
{
	int op;
	char *mne;
	r = REG_VALUE(r)*8;
	multiplier = REG_VALUE(multiplier)*8;
	if (is64_type(type))
	{
		op = SVMOPCODE_IMULQ_R_R;
		mne = "imulq";
	}
	else
	{
		op = SVMOPCODE_IMULL_R_R;
		mne = "imull";
	}
	gen1(op);
	gen1(multiplier);
	gen1(r);
	tcc_listing_code(3, -1, "%s %%r%02X, %%r%02X", mne, multiplier, r);
}

ST_FUNC void gen_div_reg(int r, int divisor, int type, int usig)
{
	int op;
	char *mne;
	r = REG_VALUE(r)*8;
	divisor = REG_VALUE(divisor)*8;
	if (is64_type(type))
	{
		if (usig)
		{
			op = SVMOPCODE_UDIVQ_R_R;
			mne = "udivq";
		}
		else
		{
			op = SVMOPCODE_DIVQ_R_R;
			mne = "divq";
		}
	}
	else
	{
		if (usig)
		{
			op = SVMOPCODE_UDIVL_R_R;
			mne = "udivl";
		}
		else
		{
			op = SVMOPCODE_DIVL_R_R;
			mne = "divl";
		}
	}
	gen1(op);
	gen1(divisor);
	gen1(r);
	tcc_listing_code(3, -1, "%s %%r%02X, %%r%02X", mne, divisor, r);
}

ST_FUNC void gen_mod_reg(int r, int divisor, int type, int usig)
{
	int op;
	char *mne;
	r = REG_VALUE(r)*8;
	divisor = REG_VALUE(divisor)*8;
	if (is64_type(type))
	{
		if (usig)
		{
			op = SVMOPCODE_UMODQ_R_R;
			mne = "umodq";
		}
		else
		{
			op = SVMOPCODE_MODQ_R_R;
			mne = "modq";
		}
	}
	else
	{
		if (usig)
		{
			op = SVMOPCODE_UMODL_R_R;
			mne = "umodl";
		}
		else
		{
			op = SVMOPCODE_MODL_R_R;
			mne = "modl";
		}
	}
	gen1(op);
	gen1(divisor);
	gen1(r);
	tcc_listing_code(3, -1, "%s %%r%02X, %%r%02X", mne, divisor, r);
}

ST_FUNC void gen_shl_reg(int r, int shift, int type)
{
	int op;
	char *mne;
	r = REG_VALUE(r)*8;
	shift = REG_VALUE(shift)*8;
	if (is64_type(type))
	{
		op = SVMOPCODE_SHLQ_R_R;
		mne = "shlq";
	}
	else
	{
		op = SVMOPCODE_SHLL_R_R;
		mne = "shll";
	}
	gen1(op);
	gen1(shift);
	gen1(r);
	tcc_listing_code(3, -1, "%s %%r%02X, %%r%02X", mne, shift, r);
}

ST_FUNC void gen_shl_const(int r, int shift, int type)
{
	int op;
	char *mne;
	r = REG_VALUE(r)*8;
	if (is64_type(type))
	{
		op = SVMOPCODE_SHLQ_IMM_R;
		mne = "shlq";
		shift &= 63;
	}
	else
	{
		op = SVMOPCODE_SHLL_IMM_R;
		mne = "shll";
		shift &= 31;
	}
	gen1(op);
	gen1(shift);
	gen1(r);
	tcc_listing_code(3, -1, "%s $%d, %%r%02X", mne, shift, r);
}

ST_FUNC void gen_shr_reg(int r, int shift, int type)
{
	int op;
	char *mne;
	r = REG_VALUE(r)*8;
	shift = REG_VALUE(shift)*8;
	if (is64_type(type))
	{
		op = SVMOPCODE_SHRQ_R_R;
		mne = "shrq";
	}
	else
	{
		op = SVMOPCODE_SHRL_R_R;
		mne = "shrl";
	}
	gen1(op);
	gen1(shift);
	gen1(r);
	tcc_listing_code(3, -1, "%s %%r%02X, %%r%02X", mne, shift, r);
}

ST_FUNC void gen_shr_const(int r, int shift, int type)
{
	int op;
	char *mne;
	r = REG_VALUE(r)*8;
	if (is64_type(type))
	{
		op = SVMOPCODE_SHRQ_IMM_R;
		mne = "shrq";
		shift &= 63;
	}
	else
	{
		op = SVMOPCODE_SHRL_IMM_R;
		mne = "shrl";
		shift &= 31;
	}
	gen1(op);
	gen1(shift);
	gen1(r);
	tcc_listing_code(3, -1, "%s $%d, %%r%02X", mne, shift, r);
}

ST_FUNC void gen_sar_reg(int r, int shift, int type)
{
	int op;
	char *mne;
	r = REG_VALUE(r)*8;
	shift = REG_VALUE(shift)*8;
	if (is64_type(type))
	{
		op = SVMOPCODE_SARQ_R_R;
		mne = "sarq";
	}
	else
	{
		op = SVMOPCODE_SARL_R_R;
		mne = "sarl";
	}
	gen1(op);
	gen1(shift);
	gen1(r);
	tcc_listing_code(3, -1, "%s %%r%02X, %%r%02X", mne, shift, r);
}

ST_FUNC void gen_sar_const(int r, int shift, int type)
{
	int op;
	char *mne;
	r = REG_VALUE(r)*8;
	if (is64_type(type))
	{
		op = SVMOPCODE_SARQ_IMM_R;
		mne = "sarq";
		shift &= 63;
	}
	else
	{
		op = SVMOPCODE_SARL_IMM_R;
		mne = "sarl";
		shift &= 31;
	}
	gen1(op);
	gen1(shift);
	gen1(r);
	tcc_listing_code(3, -1, "%s $%d, %%r%02X", mne, shift, r);
}

enum {
	ALUOP_ADD,
	ALUOP_OR,
	ALUOP_ADC,
	ALUOP_SBB,
	ALUOP_AND,
	ALUOP_SUB,
	ALUOP_XOR,
	ALUOP_CMP,
};

ST_FUNC void gen_aluop_reg(int aop, int r, int r2, int type)
{
	int op;
	char *mne;
	int ll = is64_type(type);
	r = REG_VALUE(r)*8;
	r2 = REG_VALUE(r2)*8;
	switch(aop)
	{
	case ALUOP_ADD:
		mne = ll ? "addq" : "addl";
		op = ll ? SVMOPCODE_ADDQ_R_R : SVMOPCODE_ADDL_R_R;
		break;
	case ALUOP_OR:
		mne = ll ? "orq" : "orl";
		op = ll ? SVMOPCODE_ORQ_R_R : SVMOPCODE_ORL_R_R;
		break;
	case ALUOP_ADC:
		mne = ll ? "adcq" : "adcl";
		op = ll ? SVMOPCODE_ADCQ_R_R : SVMOPCODE_ADCL_R_R;
		break;
	case ALUOP_SBB:
		mne = ll ? "sbbq" : "sbbl";
		op = ll ? SVMOPCODE_SBBQ_R_R : SVMOPCODE_SBBL_R_R;
		break;
	case ALUOP_AND:
		mne = ll ? "andq" : "andl";
		op = ll ? SVMOPCODE_ANDQ_R_R : SVMOPCODE_ANDL_R_R;
		break;
	case ALUOP_SUB:
		mne = ll ? "subq" : "subl";
		op = ll ? SVMOPCODE_SUBQ_R_R : SVMOPCODE_SUBL_R_R;
		break;
	case ALUOP_XOR:
		mne = ll ? "xorq" : "xorl";
		op = ll ? SVMOPCODE_XORQ_R_R : SVMOPCODE_XORL_R_R;
		break;
	case ALUOP_CMP:
		mne = ll ? "cmpq" : "cmpl";
		op = ll ? SVMOPCODE_CMPQ_R_R : SVMOPCODE_CMPL_R_R;
		break;
	default:
		assert(0);
	};
	gen1(op);
	gen1(r2);
	gen1(r);
	tcc_listing_code(3, -1, "%s %%r%02X, %%r%02X", mne, r2, r);
}

ST_FUNC void gen_aluop_const(int aop, int r, int64_t c, int type)
{
	int op;
	char *mne;
	int ll = is64_type(type);
	int imm8;
	ll = is64_type(type);
	if (!ll)
		c = (int32_t)c;
	imm8 = (c == (char)c);
	r = REG_VALUE(r)*8;
	switch(aop)
	{
	case ALUOP_ADD:
		mne = ll ? "addq" : "addl";
		if (imm8)
			op = ll ? SVMOPCODE_ADDQ_IM8S_R : SVMOPCODE_ADDL_IM8S_R;
		else
			op = ll ? SVMOPCODE_ADDQ_IMM_R : SVMOPCODE_ADDL_IMM_R;
		break;
	case ALUOP_OR:
		mne = ll ? "orq" : "orl";
		if (imm8)
			op = ll ? SVMOPCODE_ORQ_IM8S_R : SVMOPCODE_ORL_IM8S_R;
		else
			op = ll ? SVMOPCODE_ORQ_IMM_R : SVMOPCODE_ORL_IMM_R;
		break;
	case ALUOP_ADC:
		mne = ll ? "adcq" : "adcl";
		if (imm8)
			op = ll ? SVMOPCODE_ADCQ_IM8S_R : SVMOPCODE_ADCL_IM8S_R;
		else
			op = ll ? SVMOPCODE_ADCQ_IMM_R : SVMOPCODE_ADCL_IMM_R;
		break;
	case ALUOP_SBB:
		mne = ll ? "sbbq" : "sbbl";
		if (imm8)
			op = ll ? SVMOPCODE_SBBQ_IM8S_R : SVMOPCODE_SBBL_IM8S_R;
		else
			op = ll ? SVMOPCODE_SBBQ_IMM_R : SVMOPCODE_SBBL_IMM_R;
		break;
	case ALUOP_AND:
		mne = ll ? "andq" : "andl";
		if (imm8)
			op = ll ? SVMOPCODE_ANDQ_IM8S_R : SVMOPCODE_ANDL_IM8S_R;
		else
			op = ll ? SVMOPCODE_ANDQ_IMM_R : SVMOPCODE_ANDL_IMM_R;
		break;
	case ALUOP_SUB:
		mne = ll ? "subq" : "subl";
		if (imm8)
			op = ll ? SVMOPCODE_SUBQ_IM8S_R : SVMOPCODE_SUBL_IM8S_R;
		else
			op = ll ? SVMOPCODE_SUBQ_IMM_R : SVMOPCODE_SUBL_IMM_R;
		break;
	case ALUOP_XOR:
		mne = ll ? "xorq" : "xorl";
		if (imm8)
			op = ll ? SVMOPCODE_XORQ_IM8S_R : SVMOPCODE_XORL_IM8S_R;
		else
			op = ll ? SVMOPCODE_XORQ_IMM_R : SVMOPCODE_XORL_IMM_R;
		break;
	case ALUOP_CMP:
		mne = ll ? "cmpq" : "cmpl";
		if (imm8)
			op = ll ? SVMOPCODE_CMPQ_IM8S_R : SVMOPCODE_CMPL_IM8S_R;
		else
			op = ll ? SVMOPCODE_CMPQ_IMM_R : SVMOPCODE_CMPL_IMM_R;
		break;
	default:
		assert(0);
	};
	gen1(op);
	if (imm8)
		gen1((int)c);
	else if (ll)
		gen8(c);
	else
		gen4((int)c);
	gen1(r);
	if (imm8)
		tcc_listing_code(1 + 1 + 1, -1, "%s $%d, %%r%02X", mne, (int)c, r);
	else
		tcc_listing_code(1 + (ll ? 8 : 4) + 1, -1, "%s $%lld, %%r%02X", mne, c, r);
}

ST_FUNC void gen_fcmp(int r, int r2, int type, int unordered)
{
	int op;
	char *mne;
	int dbl = ((type & VT_BTYPE) == VT_DOUBLE || (type & VT_BTYPE) == VT_LDOUBLE);
	r = REG_VALUE(r)*8;
	r2 = REG_VALUE(r2)*8;
	if (unordered && dbl)
	{
		op = SVMOPCODE_FUCMPD_R_R;
		mne = "fucmpd";
	}
	else if (unordered && !dbl)
	{
		op = SVMOPCODE_FCMPD_R_R;
		mne = "fcmpd";
	}
	else if (!unordered && dbl)
	{
		op = SVMOPCODE_FUCMPS_R_R;
		mne = "fucmps";
	}
	else if (!unordered && !dbl)
	{
		op = SVMOPCODE_FCMPS_R_R;
		mne = "fcmps";
	}
	gen1(op);
	gen1(r2);
	gen1(r);
	tcc_listing_code(3, -1, "%s %%r%02X, %%r%02X", mne, r2, r);
}

ST_FUNC void gen_fadd(int r, int r2, int type)
{
	int op;
	char *mne;
	int dbl = ((type & VT_BTYPE) == VT_DOUBLE || (type & VT_BTYPE) == VT_LDOUBLE);
	r = REG_VALUE(r)*8;
	r2 = REG_VALUE(r2)*8;
	if (dbl)
	{
		op = SVMOPCODE_FADDD_R_R;
		mne = "faddd";
	}
	else
	{
		op = SVMOPCODE_FADDS_R_R;
		mne = "fadds";
	}
	gen1(op);
	gen1(r2);
	gen1(r);
	tcc_listing_code(3, -1, "%s %%r%02X, %%r%02X", mne, r2, r);
}

ST_FUNC void gen_fsub(int r, int r2, int type)
{
	int op;
	char *mne;
	int dbl = ((type & VT_BTYPE) == VT_DOUBLE || (type & VT_BTYPE) == VT_LDOUBLE);
	r = REG_VALUE(r)*8;
	r2 = REG_VALUE(r2)*8;
	if (dbl)
	{
		op = SVMOPCODE_FSUBD_R_R;
		mne = "fsubd";
	}
	else
	{
		op = SVMOPCODE_FSUBS_R_R;
		mne = "fsubs";
	}
	gen1(op);
	gen1(r2);
	gen1(r);
	tcc_listing_code(3, -1, "%s %%r%02X, %%r%02X", mne, r2, r);
}

ST_FUNC void gen_fmul(int r, int r2, int type)
{
	int op;
	char *mne;
	int dbl = ((type & VT_BTYPE) == VT_DOUBLE || (type & VT_BTYPE) == VT_LDOUBLE);
	r = REG_VALUE(r)*8;
	r2 = REG_VALUE(r2)*8;
	if (dbl)
	{
		op = SVMOPCODE_FMULD_R_R;
		mne = "fmuld";
	}
	else
	{
		op = SVMOPCODE_FMULS_R_R;
		mne = "fmuls";
	}
	gen1(op);
	gen1(r2);
	gen1(r);
	tcc_listing_code(3, -1, "%s %%r%02X, %%r%02X", mne, r2, r);
}

ST_FUNC void gen_fdiv(int r, int r2, int type)
{
	int op;
	char *mne;
	int dbl = ((type & VT_BTYPE) == VT_DOUBLE || (type & VT_BTYPE) == VT_LDOUBLE);
	r = REG_VALUE(r)*8;
	r2 = REG_VALUE(r2)*8;
	if (dbl)
	{
		op = SVMOPCODE_FDIVD_R_R;
		mne = "fdivd";
	}
	else
	{
		op = SVMOPCODE_FDIVS_R_R;
		mne = "fdivs";
	}
	gen1(op);
	gen1(r2);
	gen1(r);
	tcc_listing_code(3, -1, "%s %%r%02X, %%r%02X", mne, r2, r);
}

ST_FUNC void gen_cvtitof(int dst, int src, int dt, int st)
{
	int op;
	char *mne;
	int dbl = ((dt & VT_BTYPE) == VT_DOUBLE || (dt & VT_BTYPE) == VT_LDOUBLE);
	int ll = 0, ul = 0;
	dst = REG_VALUE(dst)*8;
	src = REG_VALUE(src)*8;
    if ((st & (VT_BTYPE | VT_UNSIGNED)) == (VT_INT | VT_UNSIGNED))
	{
		gen_load_reg(src | TREG_CAST64, src, st);
		ll = 1;
	}
	else if ((st & (VT_BTYPE | VT_UNSIGNED)) == (VT_LLONG | VT_UNSIGNED)) 
	{
		ul = 1;
    }
	else if ((st & VT_BTYPE) == VT_LLONG) 
	{
		ll = 1;
    }
	if (dbl)
	{
		if (ul)
		{
			op = SVMOPCODE_FCVTUQ2D_R_R;
			mne = "fcvtuq2d";
		}
		else if (ll)
		{
			op = SVMOPCODE_FCVTQ2D_R_R;
			mne = "fcvtq2d";
		}
		else
		{
			op = SVMOPCODE_FCVTL2D_R_R;
			mne = "fcvtl2d";
		}
	}
	else
	{
		if (ul)
		{
			op = SVMOPCODE_FCVTUQ2S_R_R;
			mne = "fcvtuq2s";
		}
		else if (ll)
		{
			op = SVMOPCODE_FCVTQ2S_R_R;
			mne = "fcvtq2s";
		}
		else
		{
			op = SVMOPCODE_FCVTL2S_R_R;
			mne = "fcvtl2s";
		}
	}
	gen1(op);
	gen1(src);
	gen1(dst);
	tcc_listing_code(3, -1, "%s %%r%02X, %%r%02X", mne, src, dst);
}

ST_FUNC void gen_cvtftof(int dst, int src, int dt, int st)
{
	int op;
	char *mne = NULL;
	int ll;
	dst = REG_VALUE(dst)*8;
	src = REG_VALUE(src)*8;
    if (st == VT_FLOAT)
	{
        if (dt == VT_DOUBLE || dt == VT_LDOUBLE)
		{
			op = SVMOPCODE_FCVTS2D_R_R;
			mne = "fcvts2d";
        }
    }
	else if (st == VT_DOUBLE || st == VT_LDOUBLE)
	{
        if (dt == VT_FLOAT)
		{
			op = SVMOPCODE_FCVTD2S_R_R;
			mne = "fcvtd2s";
        }
    }
	if (mne)
	{
		gen1(op);
		gen1(src);
		gen1(dst);
		tcc_listing_code(3, -1, "%s %%r%02X, %%r%02X", mne, src, dst);
	}
}

ST_FUNC void gen_cvtftoi(int dst, int src, int dt, int st)
{
	int op;
	char *mne;
	int dbl = ((st & VT_BTYPE) == VT_DOUBLE || (st & VT_BTYPE) == VT_LDOUBLE);
	int ll = 0, ul = 0;
	dst = REG_VALUE(dst)*8;
	src = REG_VALUE(src)*8;
    if ((dt & (VT_BTYPE | VT_UNSIGNED)) == (VT_INT | VT_UNSIGNED))
		ll = 1;
	else if ((dt & (VT_BTYPE | VT_UNSIGNED)) == (VT_LLONG | VT_UNSIGNED)) 
		ul = 1;
	else if ((dt & VT_BTYPE) == VT_LLONG) 
		ll = 1;
	if (dbl)
	{
		if (ul)
		{
			op = SVMOPCODE_FCVTD2UQ_R_R;
			mne = "fcvtd2uq";
		}
		else if (ll)
		{
			op = SVMOPCODE_FCVTD2Q_R_R;
			mne = "fcvtd2q";
		}
		else
		{
			op = SVMOPCODE_FCVTD2L_R_R;
			mne = "fcvtd2l";
		}
	}
	else
	{
		if (ul)
		{
			op = SVMOPCODE_FCVTS2UQ_R_R;
			mne = "fcvts2uq";
		}
		else if (ll)
		{
			op = SVMOPCODE_FCVTS2Q_R_R;
			mne = "fcvts2q";
		}
		else
		{
			op = SVMOPCODE_FCVTS2L_R_R;
			mne = "fcvts2l";
		}
	}
	gen1(op);
	gen1(src);
	gen1(dst);
	tcc_listing_code(3, -1, "%s %%r%02X, %%r%02X", mne, src, dst);
}


/* load register 'r' from value 'sv' */
void load(int r, SValue *sv)
{
	int v, t, ft, fc, fr;
	SValue v1;

	SValue v2;
	sv = pe_getimport(sv, &v2);

	fr = sv->r;
	ft = sv->type.t & ~VT_DEFSIGN;
	fc = sv->c.ul;

	v = fr & VT_VALMASK;
	if (fr & VT_LVAL)
	{
		int b, ll;
		if (v == VT_LLOCAL)
		{
			v1.type.t = VT_PTR;
			v1.r = VT_LOCAL | VT_LVAL;
			v1.c.ul = fc;
			fr = r;
			if (!(reg_classes[fr] & RC_INT))
				fr = get_reg(RC_INT);
			load(fr, &v1);
		}
		assert(((ft & VT_TYPE) == VT_BYTE || (ft & VT_TYPE) == VT_BOOL)
			|| ((ft & VT_TYPE) == (VT_BYTE | VT_UNSIGNED))
			|| ((ft & VT_TYPE) == VT_SHORT)
			|| ((ft & VT_TYPE) == (VT_SHORT | VT_UNSIGNED)) 
			|| ((ft & VT_BTYPE) == VT_INT) || ((ft & VT_BTYPE) == VT_LLONG)
			|| ((ft & VT_BTYPE) == VT_FLOAT) || ((ft & VT_BTYPE) == VT_DOUBLE)
			|| ((ft & VT_BTYPE) == VT_LDOUBLE) || ((ft & VT_BTYPE) == VT_PTR)
			|| ((ft & VT_BTYPE) == VT_ENUM)	|| ((ft & VT_BTYPE) == VT_FUNC));
		if (v == VT_CONST)
			gen_load_rip(r, fr, sv->sym, fc, ft);
		else if (v == VT_LOCAL)
			gen_load_rbp(r, fc, ft);
		else
			gen_load_ind(r, fr, ft);
	}
	else
	{
		if (v == VT_CONST)
		{
			if (fr & VT_SYM)
				gen_lea_rip(r, fr, sv->sym, fc); /* lea disp(%rip), r */
			else
				gen_load_const(r, &sv->c, ft); /* mov $xx, r */
		}
		else if (v == VT_LOCAL)
		{
			gen_lea_rbp(r, VT_LOCAL, sv->sym, fc); /* lea disp(%rbp), r */
		}
		else if (v == VT_CMP)
		{
			int t;
			if (fc & 0x100)
			{
				/* This was a float compare.  If the parity bit is
				set the result was unordered, meaning false for everything
				except TOK_NE, and true for TOK_NE.  */
				int l;
				CValue c;
				if (is64_type(ft))
					c.ull = ((fc & ~0x100) != TOK_NE) ? 0 : 1;
				else
					c.ul = ((fc & ~0x100) != TOK_NE) ? 0 : 1;
				gen_load_const(r, &c, ft); /* mov $l, r */
				fc &= ~0x100;
				t = gen_branch_psym(0, TOK_PE);
				gen_set_reg(r, fc); /* setxx %r */
				gsym(t);
			}
			else
			{
				gen_set_reg(r, fc); /* setxx %r */
			}
		}
		else if (v == VT_JMP || v == VT_JMPI)
		{
			CValue val;
			val.ul = v & 1;
			gen_load_const(r, &val, VT_INT); /* mov $1, r */
			gen_jmp_addr(ind + 2 + 3);		 /* jmp after */
			gsym(fc);
			val.ul ^= 1;
			gen_load_const(r, &val, VT_INT); /* mov $0, r */
		}
		else if (v != r)
		{
			gen_load_reg(r, v, ft);
		}
	}
}

/* store register 'r' in lvalue 'v' */
void store(int r, SValue *dv)
{
    int fr, bt, ft, fc, v;

    SValue v2;
    dv = pe_getimport(dv, &v2);

    ft = dv->type.t & ~VT_DEFSIGN;
    fc = dv->c.ul;
    fr = dv->r;
    bt = ft & VT_BTYPE;

	v = fr & VT_VALMASK;
    if (v == VT_CONST)
	{
		gen_store_rip(r, fr, dv->sym, fc, ft);
	}
	else if (v == VT_LOCAL)
	{
		gen_store_rbp(r, fc, ft);
	}
	else if (fr & VT_LVAL)
	{
		assert(REG_VALUE(fr) <= 0x1f);
		gen_store_ind(r, fr, ft);
    }
	else if (fr != r)
	{
		gen_load_reg(fr, r, ft);
    }
}

/* 'is_jmp' is '1' if it is a jump */
static void gcall_or_jmp(int is_jmp)
{
    if ((vtop->r & (VT_VALMASK | VT_LVAL)) == VT_CONST)
	{
        /* constant case */
		if (is_jmp)
			gen_jmp_rel(vtop->r, vtop->sym, vtop->c.ul);
		else
			gen_call_rel(vtop->r, vtop->sym, vtop->c.ul);
    }
	else
	{
        /* otherwise, indirect call */
        int r = gv(RC_INT);
		if (is_jmp)
			gen_jmp_reg(r);
		else
			gen_call_reg(r);
    }
}

/* Return the number of registers needed to return the struct, or 0 if
   returning via struct pointer. */
ST_FUNC int gfunc_sret(CType *vt, int variadic, CType *ret, int *ret_align)
{
    int size, align;

    *ret_align = 1; // Never have to re-align return values for SVM
    size = type_size(vt, &align);
    
	if (size > 8)
        return 0;

    ret->ref = NULL;
    ret->t = (size > 4) ? VT_LLONG : VT_INT;
    return 1;
}

static void gfunc_svmcall(int nb_args)
{
    int size, align, r, args_size, i, func_call;
    Sym *func_sym;
    
    args_size = 0;
    for(i = 0;i < nb_args; i++)
	{
        if ((vtop->type.t & VT_BTYPE) == VT_STRUCT)
		{
            size = type_size(&vtop->type, &align);
            /* align to stack align size */
            size = (size + 7) & ~7;
            /* allocate the necessary size on stack */
			gen_addsp(-size, 1);
            /* generate structure store */
            r = get_reg(RC_INT);
			gen_stsp_reg(r);
            vset(&vtop->type, r | VT_LVAL, 0);
            vswap();
            vstore();
            args_size += size;
        }
		else
		{
            /* simple type (currently always same size) */
            r = gv(is_float(vtop->type.t) ? RC_FLOAT : RC_INT);
            size = 8;
			gen_push_reg(r, vtop->type.t);
            args_size += size;
        }
        vtop--;
    }
    save_regs(0); /* save used temporary registers */
    func_sym = vtop->type.ref;
    func_call = func_sym->a.func_call;
    gcall_or_jmp(0);

    if (args_size && func_call != FUNC_STDCALL)
		gen_addsp(args_size, 1);
    vtop--;
}

static void gfunc_nativecall(int nb_args)
{
    int size, align, r, args_size, i, func_call;
    Sym *func_sym;
    
    args_size = 0;
    for(i = 0;i < nb_args; i++)
	{
        if ((vtop->type.t & VT_BTYPE) == VT_STRUCT)
		{
            size = type_size(&vtop->type, &align);
            /* align to stack align size */
            size = (size + 3) & ~3;
            /* allocate the necessary size on stack */
			gen_addsp(-size, 1);
            /* generate structure store */
            r = get_reg(RC_INT);
			gen_stsp_reg(r);
            vset(&vtop->type, r | VT_LVAL, 0);
            vswap();
            vstore();
            args_size += size;
        }
		else
		{
            /* simple type */
            r = gv(is_float(vtop->type.t) ? RC_FLOAT : RC_INT);
#ifdef TCC_TARGET_SVM64
			gen_push_reg(r, vtop->type.t);
			size = 8;
#else
            if (is64_type(vtop->type.t))
			{
				gen_push_reg(r, vtop->type.t);
				size = 8;
			}
			else
			{
				gen_push_reg32(r, vtop->type.t);
				size = 4;
			}
#endif
            args_size += size;
        }
        vtop--;
    }
    save_regs(0); /* save used temporary registers */
    func_sym = vtop->type.ref;
    func_call = func_sym->a.func_call;
	if (func_call == FUNC_NATIVETHISCALL)
		gen_ntcall_reg(gv(RC_INT), args_size);
	else
		gen_ncall_reg(gv(RC_INT), args_size);

    if (args_size)
        gen_addsp(args_size, 1);
    vtop--;
}

/* Generate function call. The function address is pushed first, then
   all the parameters in call order. This functions pops all the
   parameters and the function address. */
ST_FUNC void gfunc_call(int nb_args)
{
    int func_call = vtop[-nb_args].type.ref->a.func_call;
	if (func_call == FUNC_NATIVECALL || func_call == FUNC_NATIVETHISCALL)
		gfunc_nativecall(nb_args);
	else
		gfunc_svmcall(nb_args);
}

#define FUNC_PROLOG_SIZE 5

/* generate function prolog of type 't' */
ST_FUNC void gfunc_prolog(CType *func_type)
{
    int addr, align, size, func_call;
    int param_index, param_addr;
	int func_naked;
    Sym *sym;
    CType *type;

    sym = func_type->ref;
    func_call = sym->a.func_call;
	func_naked = sym->a.func_naked;
    addr = func_naked ? 8 : 16;
    loc = 0;
    func_vc = 0;

    param_index = 0;

	if (!func_naked)
		gen_enter(-1);

    func_sub_sp_offset = ind;
    /* if the function returns a structure, then add an
       implicit pointer parameter */
    func_vt = sym->type;
    func_var = (sym->c == FUNC_ELLIPSIS);
    size = type_size(&func_vt,&align);
    if (((func_vt.t & VT_BTYPE) == VT_STRUCT) && (size > 8))
	{
        func_vc = addr;
        addr += 8;
        param_index++;
    }
    /* define parameters */
    while ((sym = sym->next) != NULL)
	{
        type = &sym->type;
        size = type_size(type, &align);
        size = (size + 7) & ~7;
        param_addr = addr;
        addr += size;
        sym_push(sym->v & ~SYM_FIELD, type,
                 VT_LOCAL | lvalue_type(type->t), param_addr);
        param_index++;
    }
    func_ret_sub = 0;

    if (func_call == FUNC_STDCALL)
        func_ret_sub = addr - (func_naked ? 8 : 16);
}

/* generate function epilog */
ST_FUNC void gfunc_epilog(void)
{
    int v, saved_ind;

	gen_leave();
	gen_ret(func_ret_sub);

	/* align local size to qword & save local variables */    
    v = (-loc + 7) & -8; 
    saved_ind = ind;
    ind = func_sub_sp_offset - FUNC_PROLOG_SIZE;
	gen_enter(v);
    ind = saved_ind;
	tcc_listing_print("%s@__LOCAL_SIZE__=%d\n", funcname, v);
}

/* generate a jump to a label */
ST_FUNC int gjmp(int t)
{
	return gen_jmp_psym(t);
}

/* generate a jump to a fixed address */
ST_FUNC void gjmp_addr(int a)
{
	gen_jmp_addr(a);
}

/* generate a test. set 'inv' to invert test. Stack entry is popped */
ST_FUNC int gtst(int inv, int t)
{
    int v, *p;

    v = vtop->r & VT_VALMASK;
    if (v == VT_CMP)
	{
        /* fast case : can jump directly since flags are set */
		if (vtop->c.i & 0x100)
		{
			/* This was a float compare. If the parity flag is set
			the result was unordered. For anything except != this
			means false and we don't jump (anding both conditions).
			For != this means true (oring both).
			Take care about inverting the test.  We need to jump
			to our target if the result was unordered and test wasn't NE,
			otherwise if unordered we don't want to jump.  */
			vtop->c.i &= ~0x100;
			if (!inv == (vtop->c.i != TOK_NE))
			{
				int ta;
				ta = gen_branch_psym(0, TOK_PE);
				t = gen_branch_psym(t, vtop->c.i ^ inv);
				gsym(ta);
			}
			else
			{
				t = gen_branch_psym(t, TOK_PE);
				t = gen_branch_psym(t, vtop->c.i ^ inv);
			}
		}
		else
		{
			t = gen_branch_psym(t, vtop->c.i ^ inv);
		}
    }
	else
	{ 
		/* VT_JMP || VT_JMPI */
        /* && or || optimization */
        if ((v & 1) == inv)
		{
            /* insert vtop->c jump list in t */
            p = &vtop->c.i;
            while (*p != 0)
                p = (int *)(cur_text_section->data + *p);
            *p = t;
            t = vtop->c.i;
        }
		else
		{
            t = gjmp(t);
            gsym(vtop->c.i);
        }
    }
    vtop--;
    return t;
}

/* generate an integer binary operation */
void gen_opi(int op)
{
    int r, fr, ft, opc, c;
    int ll, cc;

    ll = is64_type(vtop[-1].type.t);
    cc = (vtop->r & (VT_VALMASK | VT_LVAL | VT_SYM)) == VT_CONST;

    switch(op)
	{
    case '*':
		gv2(RC_INT, RC_INT);
		gen_imul_reg(vtop[-1].r, vtop[0].r, vtop[-1].type.t);
		vtop--;
		return;

	case '/':
    case TOK_PDIV:
    case TOK_UDIV:
		gv2(RC_INT, RC_INT);
		gen_div_reg(vtop[-1].r, vtop[0].r, vtop[-1].type.t, op == TOK_UDIV);
		vtop--;
		return;

    case '%':
    case TOK_UMOD:
		gv2(RC_INT, RC_INT);
		gen_mod_reg(vtop[-1].r, vtop[0].r, vtop[-1].type.t, op == TOK_UMOD);
		vtop--;
		return;
	}

	if (cc)
	{
		vswap();
		r = gv(RC_INT);
		vswap();
		switch(op)
		{
		case TOK_SHL:
			gen_shl_const(r, vtop->c.i, vtop[-1].type.t);
			break;
		case TOK_SHR:
			gen_shr_const(r, vtop->c.i, vtop[-1].type.t);
			break;
		case TOK_SAR:
			gen_sar_const(r, vtop->c.i, vtop[-1].type.t);
			break;
		case '+':
		case TOK_ADDC1: /* add with carry generation */
			gen_aluop_const(ALUOP_ADD, r, vtop->c.ll, vtop[-1].type.t);
			break;
		case '|':
			gen_aluop_const(ALUOP_OR, r, vtop->c.ll, vtop[-1].type.t);
			break;
		case TOK_ADDC2: /* add with carry use */
			gen_aluop_const(ALUOP_ADC, r, vtop->c.ll, vtop[-1].type.t);
			break;
		case TOK_SUBC2: /* sub with carry use */
			gen_aluop_const(ALUOP_SBB, r, vtop->c.ll, vtop[-1].type.t);
			break;
		case '&':
			gen_aluop_const(ALUOP_AND, r, vtop->c.ll, vtop[-1].type.t);
			break;
		case '-':
		case TOK_SUBC1: /* sub with carry generation */
			gen_aluop_const(ALUOP_SUB, r, vtop->c.ll, vtop[-1].type.t);
			break;
	    case '^':
			gen_aluop_const(ALUOP_XOR, r, vtop->c.ll, vtop[-1].type.t);
			break;
	    default: /* compare */
			gen_aluop_const(ALUOP_CMP, r, vtop->c.ll, vtop[-1].type.t);
			break;
		}
	}
	else
	{
		gv2(RC_INT, RC_INT);
		switch(op)
		{
		case TOK_SHL:
			gen_shl_reg(vtop[-1].r, vtop[0].r, vtop[-1].type.t);
			break;
		case TOK_SHR:
			gen_shr_reg(vtop[-1].r, vtop[0].r, vtop[-1].type.t);
			break;
		case TOK_SAR:
			gen_sar_reg(vtop[-1].r, vtop[0].r, vtop[-1].type.t);
			break;
		case '+':
		case TOK_ADDC1: /* add with carry generation */
			gen_aluop_reg(ALUOP_ADD, vtop[-1].r, vtop[0].r, vtop[-1].type.t);
			break;
		case '|':
			gen_aluop_reg(ALUOP_OR, vtop[-1].r, vtop[0].r, vtop[-1].type.t);
			break;
		case TOK_ADDC2: /* add with carry use */
			gen_aluop_reg(ALUOP_ADC, vtop[-1].r, vtop[0].r, vtop[-1].type.t);
			break;
		case TOK_SUBC2: /* sub with carry use */
			gen_aluop_reg(ALUOP_SBB, vtop[-1].r, vtop[0].r, vtop[-1].type.t);
			break;
		case '&':
			gen_aluop_reg(ALUOP_AND, vtop[-1].r, vtop[0].r, vtop[-1].type.t);
			break;
		case '-':
		case TOK_SUBC1: /* sub with carry generation */
			gen_aluop_reg(ALUOP_SUB, vtop[-1].r, vtop[0].r, vtop[-1].type.t);
			break;
	    case '^':
			gen_aluop_reg(ALUOP_XOR, vtop[-1].r, vtop[0].r, vtop[-1].type.t);
			break;
	    default: /* compare */
			gen_aluop_reg(ALUOP_CMP, vtop[-1].r, vtop[0].r, vtop[-1].type.t);
			break;
		}
	}
	vtop--;
    if (op >= TOK_ULT && op <= TOK_GT)
	{
        vtop->r = VT_CMP;
        vtop->c.i = op;
    }
	return;
}

void gen_opl(int op)
{
    gen_opi(op);
}

/* generate a floating point operation 'v = t1 op t2' instruction. The
   two operands are guaranted to have the same floating point type */
void gen_opf(int op)
{
    int a, ft, fc, r;

    /* load values into registers */
    gv2(RC_FLOAT, RC_FLOAT);

    if (op >= TOK_ULT && op <= TOK_GT) 
	{
        if (op != TOK_EQ && op != TOK_NE)
		{
			if (op == TOK_LT)
                op = TOK_ULT;
			else if (op == TOK_LE)
                op = TOK_ULE;
			else if (op == TOK_GT)
                op = TOK_UGT;
			else if (op == TOK_GE)
                op = TOK_UGE;
			else
				assert(0);
        }

		gen_fcmp(vtop[-1].r, vtop[0].r, vtop->type.t, (op == TOK_EQ || op == TOK_NE));

        vtop--;
        vtop->r = VT_CMP;
		vtop->c.i = op;
		if (op == TOK_EQ || op == TOK_NE)
			vtop->c.i |= 0x100;
    }
	else
	{
        switch(op)
		{
        default:
        case '+':
			gen_fadd(vtop[-1].r, vtop[0].r, vtop->type.t);
            break;
        case '-':
			gen_fsub(vtop[-1].r, vtop[0].r, vtop->type.t);
            break;
        case '*':
			gen_fmul(vtop[-1].r, vtop[0].r, vtop->type.t);
            break;
        case '/':
			gen_fdiv(vtop[-1].r, vtop[0].r, vtop->type.t);
            break;
        }

        vtop--;
    }
}

/* convert integers to fp 't' type. Must handle 'int', 'unsigned int'
   and 'long long' cases. */
void gen_cvt_itof(int t)
{
    gv(RC_INT);
	gen_cvtitof(vtop->r, vtop->r, t, vtop->type.t);
}

/* convert from one floating point type to another */
void gen_cvt_ftof(int t)
{
    int ft, st, dt;

    st = vtop->type.t & VT_BTYPE;
    dt = t & VT_BTYPE;

	gen_cvtftof(vtop->r, vtop->r, t, vtop->type.t);
    
    if (st == VT_FLOAT)
	{
	    gv(RC_FLOAT);
        if (dt == VT_DOUBLE || dt == VT_LDOUBLE)
			gen_cvtftof(vtop->r, vtop->r, t, vtop->type.t);
    }
	else if (st == VT_DOUBLE || st == VT_LDOUBLE)
	{
        gv(RC_FLOAT);
        if (dt == VT_FLOAT)
			gen_cvtftof(vtop->r, vtop->r, t, vtop->type.t);
    }
}

/* convert fp to int 't' type */
void gen_cvt_ftoi(int t)
{
    gv(RC_FLOAT);
	gen_cvtftoi(vtop->r, vtop->r, t, vtop->type.t);
}

/* computed goto support */
ST_FUNC void ggoto(void)
{
    gcall_or_jmp(1);
    vtop--;
}

/* Save the stack pointer onto the stack */
ST_FUNC void gen_vla_sp_save(int addr)
{
	gen_stsp_rbp(addr);
}

/* Restore the SP from a location on the stack */
ST_FUNC void gen_vla_sp_restore(int addr)
{
	gen_ldsp_rbp(addr);
}

/* Subtract from the stack pointer, and push the resulting value onto the stack */
ST_FUNC void gen_vla_alloc(CType *type, int align)
{
    int r;
    r = gv(RC_INT);		/* get allocation size */    
	gen_subsp_reg(r);	/* sub %r, %rsp */
	gen_andsp_IM8S(~15);/* and $~15, %rsp - we align to 16 bytes rather than align */    
	gen_stsp_reg(r);	/* mov %rsp, %r */
    vpop();
    vset(type, r, 0);
}
