/*
 *  i386 specific functions for TCC assembler
 *
 *  Copyright (c) 2001, 2002 Fabrice Bellard
 *  Copyright (c) 2009 Frédéric Feret (x86_64 support)
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
#include "svm-opcodes.h"

#define MAX_OPERANDS 2
#define NB_SAVED_REGS 0

#define TOK_ASM_first TOK_ASM_nop
#define TOK_ASM_last TOK_ASM_vmexit
#define TOK_ASM_alllast TOK_ASM_fcvtd2s

#define OPC_JMP        0x01  /* jmp operand */
#define OPC_SHORTJMP   0x02  /* short jmp operand */
#define OPC_B          0x04  /* accepts b */
#define OPC_W          0x08  /* accepts w */
#define OPC_L          0x10  /* accepts l */
#define OPC_Q          0x20  /* accepts q */
#define OPC_DISP8S     0x40  /* -128..+127 displacement */

/* in order to compress the operand type, we use specific operands and
   we or only with EA  */
enum {
    OPT_REG,   /* warning: value is hardcoded from TOK_ASM_xxx */
    OPT_IM8,
    OPT_IM8S,
    OPT_IM16,
    OPT_IM16S,
    OPT_IM32,
    OPT_IM32S,
    OPT_IM64,
    OPT_ADDR,   /* OP_EA with only offset */
    OPT_INDIR,  /* *(expr) */
    OPT_RBP,  /* *(%rbp) */
    OPT_RSP,  /* *(%rsp) */
    OPT_RIP,  /* *(%rip) */
    /* composite types */
    OPT_COMPOSITE_FIRST,
    OPT_IM,     /* IM8 | IM8S | IM16 | IM16S | IM32 | IM64 */
    OPT_NODISP,
    /* can be ored with any OPT_xxx */
    OPT_EA = 0x80
};

#define OP_REG    (1 << OPT_REG)
#define OP_IM8    (1 << OPT_IM8)
#define OP_IM8S   (1 << OPT_IM8S)
#define OP_IM16   (1 << OPT_IM16)
#define OP_IM16S  (1 << OPT_IM16S)
#define OP_IM32   (1 << OPT_IM32)
#define OP_IM32S  (1 << OPT_IM32S)
#define OP_IM64   (1 << OPT_IM64)
#define OP_ADDR   (1 << OPT_ADDR)
#define OP_INDIR  (1 << OPT_INDIR)
#define OP_RBP    (1 << OPT_RBP)
#define OP_RSP    (1 << OPT_RSP)
#define OP_RIP    (1 << OPT_RIP)

#define OP_EA     0x40000000

typedef struct ASMInstr {
    uint16_t sym;
    uint16_t opcode;
    uint16_t instr_type;
    uint8_t nb_ops;
    uint8_t op_type[MAX_OPERANDS]; /* see OP_xxx */
} ASMInstr;

typedef struct Operand {
    uint32_t type;
    int8_t  reg; /* register, -1 if none */
    ExprValue e;
} Operand;

#define NB_TEST_OPCODES 30

static const uint8_t test_bits[NB_TEST_OPCODES] = {
 0x00, /* o */
 0x01, /* no */
 0x02, /* b */
 0x02, /* c */
 0x02, /* nae */
 0x03, /* nb */
 0x03, /* nc */
 0x03, /* ae */
 0x04, /* e */
 0x04, /* z */
 0x05, /* ne */
 0x05, /* nz */
 0x06, /* be */
 0x06, /* na */
 0x07, /* nbe */
 0x07, /* a */
 0x08, /* s */
 0x09, /* ns */
 0x0a, /* p */
 0x0a, /* pe */
 0x0b, /* np */
 0x0b, /* po */
 0x0c, /* l */
 0x0c, /* nge */
 0x0d, /* nl */
 0x0d, /* ge */
 0x0e, /* le */
 0x0e, /* ng */
 0x0f, /* nle */
 0x0f, /* g */
};

static const ASMInstr asm_instrs[] = {
#define ALT(x) x
#define DEF_ASM_OP0(name, opcode)
#define DEF_ASM_OP0L(name, opcode, instr_type) { TOK_ASM_ ## name, opcode, instr_type, 0 },
#define DEF_ASM_OP1(name, opcode, instr_type, op0) {TOK_ASM_ ## name, opcode, instr_type, 1, { op0 }},
#define DEF_ASM_OP2(name, opcode, instr_type, op0, op1) { TOK_ASM_ ## name, opcode, instr_type, 2, { op0, op1 }},
#include "svm-asm.h"
    /* last operation */
    { 0, },
};

static const uint16_t op0_codes[] = {
#define ALT(x)
#define DEF_ASM_OP0(x, opcode) opcode,
#define DEF_ASM_OP0L(name, opcode, instr_type)
#define DEF_ASM_OP1(name, opcode, instr_type, op0)
#define DEF_ASM_OP2(name, opcode, instr_type, op0, op1)
#include "svm-asm.h"
};

static int asm_parse_reg(void)
{
    int reg = 0;
    if (tok != '%')
        goto error_reg;
    next();
    if (tok >= TOK_ASM_r00 && tok <= TOK_ASM_rip)
	{
        reg = tok - TOK_ASM_r00;
    }
	else
	{
    error_reg:
        expect("register");
    }
    next();
    return reg;
}

static void parse_operand(TCCState *s1, Operand *op)
{
    ExprValue e;
    int reg, indir;
    const char *p;
	int disp = 0;

    indir = 0;
    if (tok == '*')
	{
        next();
        indir = OP_INDIR;
    }

    if (tok == '%')
	{
        next();
        if (tok >= TOK_ASM_r00 && tok <= TOK_ASM_rF8)
		{
            op->type = OP_REG;
            op->reg = tok - TOK_ASM_r00;
        }
        else if (tok == TOK_ASM_rsp)
		{
            op->type = OP_RSP;
            op->reg = tok - TOK_ASM_r00;
		}
		else
		{
        reg_error:
            tcc_error("unknown register");
        }
        next();
    } 
	else if (tok == '$')
	{
        /* constant value */
        next();
        asm_expr(s1, &e);
        op->e.v = e.v;
        op->e.sym = e.sym;
        if (!op->e.sym)
		{
            if (op->e.v == (int8_t)op->e.v)
                op->type = OP_IM8S | OP_IM8 | OP_IM16S | OP_IM16 | OP_IM32S | OP_IM32 | OP_IM64;
            else if (op->e.v == (uint8_t)op->e.v)
                op->type = OP_IM8 | OP_IM16S | OP_IM16 | OP_IM32S | OP_IM32 | OP_IM64;
            else if (op->e.v == (int16_t)op->e.v)
                op->type = OP_IM16S | OP_IM16 | OP_IM32S | OP_IM32 | OP_IM64;
            else if (op->e.v == (uint16_t)op->e.v)
                op->type = OP_IM16 | OP_IM32S | OP_IM32 | OP_IM64;
            else if (op->e.v == (int32_t)op->e.v)
                op->type = OP_IM32S | OP_IM32 | OP_IM64;
            else if (op->e.v == (uint32_t)op->e.v)
                op->type = OP_IM32 | OP_IM64;
            else if (op->e.v == (uint64_t)op->e.v)
                op->type = OP_IM64;
        }
		else
		{
	        op->type = (PTR_SIZE == 8) ? OP_IM64 : OP_IM32;
		}
    }
	else
	{
        /* address(reg) with all variants */
        op->type = OP_EA;
        op->reg = -1;
        if (tok != '(') 
		{
            asm_expr(s1, &e);
            op->e.v = e.v;
            op->e.sym = e.sym;
			disp = 1;
        }
		else
		{
            next();
            if (tok == '%')
			{
                unget_tok('(');
                op->e.v = 0;
                op->e.sym = NULL;
            }
			else
			{
                /* bracketed offset expression */
                asm_expr(s1, &e);
                if (tok != ')')
                    expect(")");
                next();
                op->e.v = e.v;
                op->e.sym = e.sym;
            }
        }
        if (tok == '(')
		{
            next();
			if (tok != '%')
				goto error_reg;
			next();
			if (tok >= TOK_ASM_r00 && tok <= TOK_ASM_rF8)
			{
				op->type |= OP_REG;
				op->reg = tok - TOK_ASM_r00;
			}
			else if (tok == TOK_ASM_rbp)
			{
				op->type |= OP_RBP;
				op->reg = tok - TOK_ASM_r00;
			}
			else if (tok == TOK_ASM_rip)
			{
				op->type |= OP_RIP;
				op->reg = tok - TOK_ASM_r00;
			}
			else
			{
			error_reg:
				expect("register");
			}
			next();
            skip(')');
        }
		if (disp && (op->type & OP_REG))
			expect("no displacement");
        if (op->reg == -1)
            op->type = OP_ADDR;
    }
    op->type |= indir;
}

/* XXX: unify with C code output ? */
ST_FUNC void gen_expr32(ExprValue *pe)
{
    gen_addr32(pe->sym ? VT_SYM : 0, pe->sym, pe->v);
}

static void gen_expr64(ExprValue *pe)
{
    gen_addr64(pe->sym ? VT_SYM : 0, pe->sym, pe->v);
}

/* XXX: unify with C code output ? */
static void gen_disp32(ExprValue *pe)
{
    Sym *sym = pe->sym;
    if (sym && sym->r == cur_text_section->sh_num)
	{
        /* same section: we can output an absolute value. Note
           that the TCC compiler behaves differently here because
           it always outputs a relocation to ease (future) code
           elimination in the linker */
        gen4(pe->v + sym->jnext - ind - 4);
    }
	else
	{
        if (sym && sym->type.t == VT_VOID)
		{
            sym->type.t = VT_FUNC;
            sym->type.ref = NULL;
        }
        gen_addrpc32(VT_SYM, sym, pe->v);
    }
}

ST_FUNC void asm_opcode(TCCState *s1, int opcode)
{
    const ASMInstr *pa;
    int i, modrm_index, reg, v, op1, is_short_jmp, disp8s;
    int nb_ops, s;
    Operand ops[MAX_OPERANDS], *pop;
    int op_type[3]; /* decoded op type */
	const char *mne;
	CString instr;
	int ind1;
	int symoffset = -1;
	char buf[256];

	ind1 = ind;

    /* get operands */
    pop = ops;
    nb_ops = 0;
    for(;;)
	{
        if (tok == ';' || tok == TOK_LINEFEED)
            break;
        if (nb_ops >= MAX_OPERANDS)
            tcc_error("incorrect number of operands");
        parse_operand(s1, pop);
        pop++;
        nb_ops++;
		if (tok == ';' || tok == TOK_LINEFEED)
			break;
        if (tok != ',')
            tcc_error("bad operand");
        next();
    }

    is_short_jmp = 0;
	disp8s = 0;
    s = 0; /* avoid warning */

    /* optimize matching by using a lookup table (no hashing is needed
       !) */
    for(pa = asm_instrs; pa->sym != 0; pa++)
	{
        if (pa->sym != opcode)
            continue;
		if (pa->instr_type & OPC_B)
			s = 1;
		else if (pa->instr_type & OPC_W)
			s = 2;
		else if (pa->instr_type & OPC_L)
			s = 3;
		else if (pa->instr_type & OPC_Q)
			s = 4;
		else
			s = 3;
        if (pa->nb_ops != nb_ops)
            continue;
        /* now decode and check each operand */
        for(i = 0; i < nb_ops; i++)
		{
            int op1, op2;
            op1 = pa->op_type[i];
            op2 = op1 & 0x1f;
            switch(op2) {
            case OPT_IM:
                v = OP_IM8 | OP_IM16 | OP_IM32 | OP_IM64;
                break;
            default:
                v = 1 << op2;
                break;
            }
            if (v == OP_INDIR && ops[i].reg != -1)
                goto next;
            if (op1 & OPT_EA)
			{
				if ((pa->instr_type & OPC_DISP8S) && ops[i].e.v != (char)ops[i].e.v)
					goto next;
				if ((ops[i].type & v) == 0)
					goto next;
				if ((ops[i].type & OP_EA) == 0)
					goto next;
                v |= OP_EA;
			}
			else
			{
				if ((ops[i].type & OP_EA) != 0)
					goto next;
		        if ((ops[i].type & v) == 0)
					goto next;
			}
            op_type[i] = v;
        }
        /* all is matching ! */
        break;
    next: ;
    }
    if (pa->sym == 0) {
        if (opcode >= TOK_ASM_first && opcode <= TOK_ASM_last)
		{
            int b;
            b = op0_codes[opcode - TOK_ASM_first];
            gen1(b);
			tcc_listing_code(ind - ind1, -1, "%s", get_tok_str(opcode, NULL));
            return;
        }
		else if (opcode <= TOK_ASM_alllast)
		{
            tcc_error("bad operand with opcode '%s'",
                  get_tok_str(opcode, NULL));
        }
		else
		{
            tcc_error("unknown opcode '%s'",
                  get_tok_str(opcode, NULL));
        }
    }

    /* now generates the operation */
    v = pa->opcode;
    if (pa->instr_type & OPC_SHORTJMP)
	{
        Sym *sym;
        int jmp_disp;

        /* see if we can really generate the jump with a byte offset */
        sym = ops[0].e.sym;
        if (!sym)
            goto no_short_jump;
        if (sym->r != cur_text_section->sh_num)
            goto no_short_jump;
        jmp_disp = ops[0].e.v + sym->jnext - ind - 2;
        if (jmp_disp == (int8_t)jmp_disp)
		{
            /* OK to generate jump */
            is_short_jmp = 1;
            ops[0].e.v = jmp_disp;
        }
		else
		{
        no_short_jump:
            if (pa->instr_type & OPC_JMP)
			{
                /* long jump will be allowed. need to modify the
                   opcode slightly */
                if (v == SVMOPCODE_JMPS)
                    v = SVMOPCODE_JMP;
            }
			else
			{
                tcc_error("invalid displacement");
            }
        }
    }
    gen1(v);
	cstr_new(&instr);
	cstr_cat(&instr, get_tok_str(opcode, NULL));
    
    /* emit operands */
    for(i = 0;i < nb_ops; i++)
	{
		cstr_cat(&instr, i ? ", " : " ");
        v = op_type[i];
        if (v & (OP_IM8 | OP_IM16 | OP_IM32 | OP_IM64 | OP_IM8S | OP_IM16S | OP_IM32S | OP_ADDR))
		{
            /* if multiple sizes are given it means we must look
               at the op size */
            if ((v | OP_IM8 | OP_IM64) == (OP_IM8 | OP_IM16 | OP_IM32 | OP_IM64))
			{
                if (s == 1)
                    v = ops[i].type & (OP_IM8 | OP_IM8S);
                else if (s == 2)
                    v = ops[i].type & (OP_IM16 | OP_IM16S);
                else if (s == 3 || (v & OP_IM64) == 0)
                    v = ops[i].type & (OP_IM32 | OP_IM32S);
                else
                    v = OP_IM64;
            }
			if (ops[i].e.sym)
			{
				symoffset = is_short_jmp ? -1 : ind - ind1;
				cstr_cat(&instr, get_tok_str(ops[i].e.sym->v, NULL));
			}
			else
			{
	            if (v & OP_IM8S)
					sprintf(buf, "$%d", (int)(char)ops[i].e.v);
	            else if (v & OP_IM16S)
					sprintf(buf, "$%d", (int)(short)ops[i].e.v);
	            else if (v & OP_IM32S)
					sprintf(buf, "$%d", (int)ops[i].e.v);
	            else if (v & OP_IM8)
					sprintf(buf, "$%u", (unsigned int)(unsigned char)ops[i].e.v);
	            else if (v & OP_IM16)
					sprintf(buf, "$%u", (unsigned int)(unsigned short)ops[i].e.v);
	            else if (v & OP_IM32)
					sprintf(buf, "$%u", (unsigned int)ops[i].e.v);
				else if (v & OP_IM64)
					sprintf(buf, "$%lldLL", (int64_t)ops[i].e.v);
				else
					sprintf(buf, "$%d", (int)ops[i].e.v);
				cstr_cat(&instr, buf);
			}
			if (v & (OP_IM8 | OP_IM8S))
			{
                if (ops[i].e.sym)
                    tcc_error("cannot relocate");
                gen1(ops[i].e.v);
            }
			else if (v & (OP_IM16 | OP_IM16S))
			{
                if (ops[i].e.sym)
                    tcc_error("cannot relocate");
                else
                    gen2(ops[i].e.v);
            }
			else
			{
                if (pa->instr_type & (OPC_JMP | OPC_SHORTJMP))
				{
                    if (is_short_jmp)
                        gen1(ops[i].e.v);
                    else
                        gen_disp32(&ops[i].e);
                }
				else
				{
					if (v & OP_ADDR)
					{
						if (PTR_SIZE == 8)
							gen_expr64(&ops[i].e);
						else
							gen_expr32(&ops[i].e);
					}
					else 
					{
						if (pa->instr_type & OPC_Q)
						{
							if ((PTR_SIZE != 8) && ops[i].e.sym)
								tcc_error("cannot relocate");
							gen_expr64(&ops[i].e);
						}
						else
						{
							if ((PTR_SIZE == 8) && ops[i].e.sym)
								tcc_error("cannot relocate");
							gen_expr32(&ops[i].e);
						}
					}
                }
            }
        }
        else if (v & OP_INDIR)
		{
			if (ops[i].reg != -1)
			{
				gen1(ops[i].reg*8);
				sprintf(buf, "(%%r%02X)", ops[i].reg*8);
				cstr_cat(&instr, buf);
			}
			else
			{
				cstr_ccat(&instr, '*');
				if (ops[i].e.sym)
				{
					symoffset = ind - ind1;
					cstr_cat(&instr, get_tok_str(ops[i].e.sym->v, NULL));
				}
				else
				{
					sprintf(buf, "$%u", (unsigned int)ops[i].e.v);
					cstr_cat(&instr, buf);
				}
				gen_disp32(&ops[i].e);
			}
		}
        else if (v & OP_EA)
		{
			if (ops[i].type & OP_REG)
			{
				gen1(ops[i].reg*8);
				sprintf(buf, "(%%r%02X)", ops[i].reg*8);
				cstr_cat(&instr, buf);
			}
			else
			{
				if (v & OP_RBP)
				{
					if (ops[i].e.sym)
						tcc_error("cannot relocate");
					if (ops[i].e.v)
					{
						sprintf(buf, "%d", (int)ops[i].e.v);
						cstr_cat(&instr, buf);
					}
					cstr_cat(&instr, "(%rbp)");
					if (pa->instr_type & OPC_DISP8S)
						gen1(ops[i].e.v);
					else
						gen_expr32(&ops[i].e);
				}
				else if (v & OP_RIP)
				{
					if (ops[i].e.sym)
					{
						symoffset = ind - ind1;
						cstr_cat(&instr, get_tok_str(ops[i].e.sym->v, NULL));
					}
					else if (ops[i].e.v)
					{
						sprintf(buf, "%d", (int)ops[i].e.v);
						cstr_cat(&instr, buf);
					}
					cstr_cat(&instr, "(%rip)");	
					gen_disp32(&ops[i].e);
				}
				else
				{
					if (ops[i].e.sym)
					{
						symoffset = ind - ind1;
						cstr_cat(&instr, get_tok_str(ops[i].e.sym->v, NULL));
					}
					else
					{
						sprintf(buf, "%d", (int)ops[i].e.v);
						cstr_cat(&instr, buf);
					}
					if ((PTR_SIZE == 8) && ops[i].e.sym)
						tcc_error("cannot relocate");
					gen_expr32(&ops[i].e);
				}
			}
		}
		else if (v & OP_RSP)
		{
			cstr_cat(&instr, "%rsp");
		}
		else if (v & OP_REG)
		{
            gen1(ops[i].reg*8);
			sprintf(buf, "%%r%02X", ops[i].reg*8);
			cstr_cat(&instr, buf);
        }
    }

	cstr_ccat(&instr, '\0');
	tcc_listing_code(ind - ind1, symoffset, "%s", instr.data);
	cstr_free(&instr);
}

/* return the constraint priority (we allocate first the lowest
   numbered constraints) */
static inline int constraint_priority(const char *str)
{
    int priority, c, pr;

    /* we take the lowest priority */
    priority = 0;
    for(;;) {
        c = *str;
        if (c == '\0')
            break;
        str++;
        switch(c) {
        case 'r':
            pr = 3;
            break;
        case 'M':
        case 'N':
        case 'I':
        case 'i':
        case 'm':
        case 'g':
            pr = 4;
            break;
        default:
            tcc_error("unknown constraint '%c'", c);
            pr = 0;
        }
        if (pr > priority)
            priority = pr;
    }
    return priority;
}

static const char *skip_constraint_modifiers(const char *p)
{
    while (*p == '=' || *p == '&' || *p == '+' || *p == '%')
        p++;
    return p;
}

#define REG_OUT_MASK 0x01
#define REG_IN_MASK  0x02

#define is_reg_allocated(reg) (regs_allocated[reg] & reg_mask)

ST_FUNC void asm_compute_constraints(ASMOperand *operands,
                                    int nb_operands, int nb_outputs,
                                    const uint8_t *clobber_regs,
                                    int *pout_reg)
{
    ASMOperand *op;
    int sorted_op[MAX_ASM_OPERANDS];
    int i, j, k, p1, p2, tmp, reg, c, reg_mask;
    const char *str;
    uint8_t regs_allocated[NB_ASM_REGS];

    /* init fields */
    for(i=0;i<nb_operands;i++)
	{
        op = &operands[i];
        op->input_index = -1;
        op->ref_index = -1;
        op->reg = -1;
        op->is_memory = 0;
        op->is_rw = 0;
    }
    /* compute constraint priority and evaluate references to output
       constraints if input constraints */
    for(i=0;i<nb_operands;i++)
	{
        op = &operands[i];
        str = op->constraint;
        str = skip_constraint_modifiers(str);
        if (isnum(*str) || *str == '[')
		{
            /* this is a reference to another constraint */
            k = find_constraint(operands, nb_operands, str, NULL);
            if ((unsigned)k >= i || i < nb_outputs)
                tcc_error("invalid reference in constraint %d ('%s')",
                      i, str);
            op->ref_index = k;
            if (operands[k].input_index >= 0)
                tcc_error("cannot reference twice the same operand");
            operands[k].input_index = i;
            op->priority = 5;
        }
		else
		{
            op->priority = constraint_priority(str);
        }
    }

    /* sort operands according to their priority */
    for(i=0;i<nb_operands;i++)
        sorted_op[i] = i;
    for(i=0;i<nb_operands - 1;i++)
	{
        for(j=i+1;j<nb_operands;j++)
		{
            p1 = operands[sorted_op[i]].priority;
            p2 = operands[sorted_op[j]].priority;
            if (p2 < p1)
			{
                tmp = sorted_op[i];
                sorted_op[i] = sorted_op[j];
                sorted_op[j] = tmp;
            }
        }
    }

    for(i = 0;i < NB_ASM_REGS; i++)
	{
        if (clobber_regs[i])
            regs_allocated[i] = REG_IN_MASK | REG_OUT_MASK;
        else
            regs_allocated[i] = 0;
    }
    /* esp cannot be used */
    regs_allocated[4] = REG_IN_MASK | REG_OUT_MASK;
    /* ebp cannot be used yet */
    regs_allocated[5] = REG_IN_MASK | REG_OUT_MASK;

    /* allocate registers and generate corresponding asm moves */
    for(i=0;i<nb_operands;i++)
	{
        j = sorted_op[i];
        op = &operands[j];
        str = op->constraint;
        /* no need to allocate references */
        if (op->ref_index >= 0)
            continue;
        /* select if register is used for output, input or both */
        if (op->input_index >= 0)
		{
            reg_mask = REG_IN_MASK | REG_OUT_MASK;
        }
		else if (j < nb_outputs)
		{
            reg_mask = REG_OUT_MASK;
        }
		else
		{
            reg_mask = REG_IN_MASK;
        }
    try_next:
        c = *str++;
        switch(c)
		{
        case '=':
            goto try_next;
        case '+':
            op->is_rw = 1;
            /* FALL THRU */
        case '&':
            if (j >= nb_outputs)
                tcc_error("'%c' modifier can only be applied to outputs", c);
            reg_mask = REG_IN_MASK | REG_OUT_MASK;
            goto try_next;
        alloc_reg:
            if (is_reg_allocated(reg))
                goto try_next;
            goto reg_found;
        case 'r':
            /* any general register */
            for(reg = 0; reg < NB_ASM_REGS; reg++) {
                if (!is_reg_allocated(reg))
                    goto reg_found;
            }
            goto try_next;
        reg_found:
            /* now we can reload in the register */
            op->reg = reg;
            regs_allocated[reg] |= reg_mask;
            break;
        case 'i':
            if (!((op->vt->r & (VT_VALMASK | VT_LVAL)) == VT_CONST))
                goto try_next;
            break;
        case 'I':
        case 'N':
        case 'M':
            if (!((op->vt->r & (VT_VALMASK | VT_LVAL | VT_SYM)) == VT_CONST))
                goto try_next;
            break;
        case 'm':
        case 'g':
            /* nothing special to do because the operand is already in
               memory, except if the pointer itself is stored in a
               memory variable (VT_LLOCAL case) */
            /* XXX: fix constant case */
            /* if it is a reference to a memory zone, it must lie
               in a register, so we reserve the register in the
               input registers and a load will be generated
               later */
            if (j < nb_outputs || c == 'm')
			{
                if ((op->vt->r & VT_VALMASK) == VT_LLOCAL)
				{
                    /* any general register */
                    for(reg = 0; reg < 8; reg++)
					{
                        if (!(regs_allocated[reg] & REG_IN_MASK))
                            goto reg_found1;
                    }
                    goto try_next;
                reg_found1:
                    /* now we can reload in the register */
                    regs_allocated[reg] |= REG_IN_MASK;
                    op->reg = reg;
                    op->is_memory = 1;
                }
            }
            break;
        default:
            tcc_error("asm constraint %d ('%s') could not be satisfied",
                  j, op->constraint);
            break;
        }
        /* if a reference is present for that operand, we assign it too */
        if (op->input_index >= 0)
            operands[op->input_index].reg = op->reg;
    }

    /* compute out_reg. It is used to store outputs registers to memory
       locations references by pointers (VT_LLOCAL case) */
    *pout_reg = -1;
    for(i=0;i<nb_operands;i++)
	{
        op = &operands[i];
        if (op->reg >= 0 &&
            (op->vt->r & VT_VALMASK) == VT_LLOCAL  &&
            !op->is_memory)
		{
            for(reg = 0; reg < 8; reg++)
			{
                if (!(regs_allocated[reg] & REG_OUT_MASK))
                    goto reg_found2;
            }
            tcc_error("could not find free output register for reloading");
        reg_found2:
            *pout_reg = reg;
            break;
        }
    }

    /* print sorted constraints */
#ifdef ASM_DEBUG
    for(i=0;i<nb_operands;i++)
	{
        j = sorted_op[i];
        op = &operands[j];
        printf("%%%d [%s]: \"%s\" r=0x%04x reg=%d\n",
               j,
               op->id ? get_tok_str(op->id, NULL) : "",
               op->constraint,
               op->vt->r,
               op->reg);
    }
    if (*pout_reg >= 0)
        printf("out_reg=%d\n", *pout_reg);
#endif
}

ST_FUNC void subst_asm_operand(CString *add_str,
                              SValue *sv, int modifier)
{
    int r, reg, size, val;
    char buf[64];

    r = sv->r;
    if ((r & VT_VALMASK) == VT_CONST)
	{
        if (!(r & VT_LVAL) && modifier != 'c' && modifier != 'n')
            cstr_ccat(add_str, '$');
        if (r & VT_SYM)
		{
            cstr_cat(add_str, get_tok_str(sv->sym->v, NULL));
            if (sv->c.i != 0)
			{
                cstr_ccat(add_str, '+');
            }
			else
			{
                return;
            }
        }
        val = sv->c.i;
        if (modifier == 'n')
            val = -val;
        snprintf(buf, sizeof(buf), "%d", sv->c.i);
        cstr_cat(add_str, buf);
    }
	else if ((r & VT_VALMASK) == VT_LOCAL)
	{
        snprintf(buf, sizeof(buf), "%d(%%rbp)", sv->c.i);
        cstr_cat(add_str, buf);
    }
	else if (r & VT_LVAL)
	{
        reg = r & VT_VALMASK;
        if (reg >= VT_CONST)
            tcc_error("internal compiler error");
        snprintf(buf, sizeof(buf), "(%%%s)",get_tok_str(TOK_ASM_r00 + reg, NULL));
        cstr_cat(add_str, buf);
    }
	else
	{
        /* register case */
        reg = r & VT_VALMASK;
        if (reg >= VT_CONST)
            tcc_error("internal compiler error");

        reg = TOK_ASM_r00 + reg;
        snprintf(buf, sizeof(buf), "%%%s", get_tok_str(reg, NULL));
        cstr_cat(add_str, buf);
    }
}

/* generate prolog and epilog code for asm statement */
ST_FUNC void asm_gen_code(ASMOperand *operands, int nb_operands,
                         int nb_outputs, int is_output,
                         uint8_t *clobber_regs,
                         int out_reg)
{
    uint8_t regs_allocated[NB_ASM_REGS];
    ASMOperand *op;
    int i, reg;
#if NB_SAVED_REGS != 0
    static uint8_t reg_saved[NB_SAVED_REGS] = { 0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31 };
#endif

    /* mark all used registers */
    memcpy(regs_allocated, clobber_regs, sizeof(regs_allocated));
    for(i = 0; i < nb_operands;i++)
	{
        op = &operands[i];
        if (op->reg >= 0)
            regs_allocated[op->reg] = 1;
    }
    if (!is_output)
	{
#if NB_SAVED_REGS != 0
        /* generate reg save code */
        for(i = 0; i < NB_SAVED_REGS; i++)
		{
            reg = reg_saved[i];
            if (regs_allocated[reg])
			{
                g(0x50 + reg);
            }
        }
#endif

        /* generate load code */
        for(i = 0; i < nb_operands; i++)
		{
            op = &operands[i];
            if (op->reg >= 0)
			{
                if ((op->vt->r & VT_VALMASK) == VT_LLOCAL &&
                    op->is_memory)
				{
                    /* memory reference case (for both input and
                       output cases) */
                    SValue sv;
                    sv = *op->vt;
                    sv.r = (sv.r & ~VT_VALMASK) | VT_LOCAL;
                    load(op->reg, &sv);
                }
				else if (i >= nb_outputs || op->is_rw)
				{
                    /* load value in register */
                    load(op->reg, op->vt);
                }
            }
        }
    } else {
        /* generate save code */
        for(i = 0 ; i < nb_outputs; i++)
		{
            op = &operands[i];
            if (op->reg >= 0)
			{
                if ((op->vt->r & VT_VALMASK) == VT_LLOCAL)
				{
                    if (!op->is_memory)
					{
                        SValue sv;
                        sv = *op->vt;
                        sv.r = (sv.r & ~VT_VALMASK) | VT_LOCAL;
                        load(out_reg, &sv);

                        sv.r = (sv.r & ~VT_VALMASK) | out_reg;
                        store(op->reg, &sv);
                    }
                }
				else
				{
                    store(op->reg, op->vt);
                }
            }
        }
#if NB_SAVED_REGS != 0
        /* generate reg restore code */
        for(i = NB_SAVED_REGS - 1; i >= 0; i--)
		{
            reg = reg_saved[i];
            if (regs_allocated[reg])
			{
                g(0x58 + reg);
            }
        }
#endif
    }
}

ST_FUNC void asm_clobber(uint8_t *clobber_regs, const char *str)
{
    int reg;
    TokenSym *ts;

    if (!strcmp(str, "memory") ||
        !strcmp(str, "cc"))
        return;
    ts = tok_alloc(str, strlen(str));
    reg = ts->tok;
    if (reg >= TOK_ASM_r00 && reg <= TOK_ASM_rF8)
	{
        reg -= TOK_ASM_r00;
    }
	else
	{
        tcc_error("invalid clobber register '%s'", str);
    }
    clobber_regs[reg] = 1;
}
