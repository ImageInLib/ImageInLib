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


/* number of available registers */
#define NB_REGS         32
#define NB_ASM_REGS     32

#define REG_VALUE(reg) ((reg) & 0x1f)

/* a register can belong to several classes. The classes must be
   sorted from more general to more precise (see gv2() code which does
   assumptions on it). */
#define RC_INT     0x0001 /* generic integer register */
#define RC_FLOAT   0x0002 /* generic float register */
#define RC_R00     0x0004
#define RC_R08     0x0008
#define RC_IRET    RC_R00 /* function return: integer register */
#define RC_LRET    RC_R08 /* function return: second integer register */
#define RC_FRET    RC_R00 /* function return: float register */
#define RC_QRET    RC_R08 /* function return: second float register */

/* pretty names for the registers */
enum {
    TREG_R00 = 0,
    TREG_R08 = 1,
    TREG_R10 = 2,
    TREG_R18 = 3,
    TREG_R20 = 4,
    TREG_R28 = 5,
    TREG_R30 = 6,
    TREG_R38 = 7,
    TREG_R40 = 8,
    TREG_R48 = 9,
    TREG_R50 = 10,
    TREG_R58 = 11,
    TREG_R60 = 12,
    TREG_R68 = 13,
    TREG_R70 = 14,
    TREG_R78 = 15,
    TREG_R80 = 16,
    TREG_R88 = 17,
    TREG_R90 = 18,
    TREG_R98 = 19,
    TREG_RA0 = 20,
    TREG_RA8 = 21,
    TREG_RB0 = 22,
    TREG_RB8 = 23,
    TREG_RC0 = 24,
    TREG_RC8 = 25,
    TREG_RD0 = 26,
    TREG_RD8 = 27,
    TREG_RE0 = 28,
    TREG_RE8 = 29,
    TREG_RF0 = 30,
    TREG_RF8 = 31,
    TREG_CAST64 = 0x20,
};

/* return registers for function */
#define REG_IRET TREG_R00 /* single word int return register */
#define REG_LRET TREG_R08 /* second word return register */
#define REG_FRET TREG_R00 /* float return register */
#define REG_QRET TREG_R08 /* second float return register */

/* pointer size, in bytes */
extern const int PTR_SIZE;

/* long double size and alignment, in bytes */
#define LDOUBLE_SIZE  8
#define LDOUBLE_ALIGN 8
/* maximum alignment (for aligned attribute support) */
#define MAX_ALIGN     8


/******************************************************/
/* ELF defines */

extern const int EM_TCC_TARGET;

/* relocation type for 32 bit data relocation */
extern const int R_PC_32;
extern const int R_DATA_32;
extern const int R_DATA_PTR;
extern const int R_JMP_SLOT;
extern const int R_COPY;

extern const int ELF_START_ADDR;
extern const int ELF_PAGE_SIZE;

