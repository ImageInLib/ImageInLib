// movl %s, %d					- 00 <s> <d>
#define SVMOPCODE_MOVL_R_R		0x00
// movl $imm32, %r				- 01 <i> <i> <i> <i> <r>
#define SVMOPCODE_MOVL_IMM_R	0x01
// movl disp(%rip), r			- 02 <d> <d> <d> <d> <r>
#define SVMOPCODE_MOVL_RIP_R	0x02
// movl disp8(%rbp), r			- 03 <d> <r>
#define SVMOPCODE_MOVL_RBP8_R	0x03
// movl disp(%rbp), r			- 04 <d> <d> <d> <d> <r>
#define SVMOPCODE_MOVL_RBP_R	0x04
// movl (%s), %d				- 05 <s> <d>
#define SVMOPCODE_MOVL_IND_R	0x05

// movl $imm8, %r				- 06 <i> <r>
#define SVMOPCODE_MOVL_IM8S_R	0x06
// movl $imm16, %r				- 07 <i> <i> <r>
#define SVMOPCODE_MOVL_IM16S_R	0x07

// movq %s, %d					- 08 <s> <d>
#define SVMOPCODE_MOVQ_R_R		0x08
// movq $imm64, %r				- 09 <i> <i> <i> <i> <i> <i> <i> <i> <r>
#define SVMOPCODE_MOVQ_IMM_R	0x09
// movq disp(%rip), r			- 0a <d> <d> <d> <d> <r>
#define SVMOPCODE_MOVQ_RIP_R	0x0a
// movq disp8(%rbp), r			- 0b <d> <r>
#define SVMOPCODE_MOVQ_RBP8_R	0x0b
// movq disp(%rbp), r			- 0c <d> <d> <d> <d> <r>
#define SVMOPCODE_MOVQ_RBP_R	0x0c
// movq (%s), %d				- 0d <s> <d>
#define SVMOPCODE_MOVQ_IND_R	0x0d

// call %r						- 0e <r>
#define SVMOPCODE_CALL_R		0x0e
// call $rel					- 0f <i> <i> <i> <i>
#define SVMOPCODE_CALL			0x0f

// movb r, disp(%rip)			- 10 <r> <d> <d> <d> <d>
#define SVMOPCODE_MOVB_R_RIP	0x10
// movb r, disp8(%rbp)			- 11 <r> <d>
#define SVMOPCODE_MOVB_R_RBP8	0x11
// movb r, disp(%rbp)			- 12 <r> <d> <d> <d> <d>
#define SVMOPCODE_MOVB_R_RBP	0x12
// movb %s, (%d)				- 13 <s> <d>
#define SVMOPCODE_MOVB_R_IND	0x13
// movsbl %s, %d				- 14 <s> <d>
#define SVMOPCODE_MOVSBL_R_R	0x14
// movsbl disp(%rip), %r		- 15 <d> <d> <d> <d> <r>
#define SVMOPCODE_MOVSBL_RIP_R	0x15
// movsbl disp(%rbp), r			- 16 <d> <d> <d> <d> <r>
#define SVMOPCODE_MOVSBL_RBP_R	0x16
// movsbl (%s), %d				- 17 <s> <d>
#define SVMOPCODE_MOVSBL_IND_R	0x17
// movzbl %s, %d				- 18 <s> <d>
#define SVMOPCODE_MOVZBL_R_R	0x18
// movzbl disp(%rip), %r		- 19 <d> <d> <d> <d> <r>
#define SVMOPCODE_MOVZBL_RIP_R	0x19
// movzbl disp(%rbp), r			- 1a <d> <d> <d> <d> <r>
#define SVMOPCODE_MOVZBL_RBP_R	0x1a
// movzbl (%s), %d				- 1b <s> <d>
#define SVMOPCODE_MOVZBL_IND_R	0x1b

// pushl %r						- 1c <r>
#define SVMOPCODE_PUSHL_R		0x1c
// pushq %r						- 1d <r>
#define SVMOPCODE_PUSHQ_R		0x1d
// push32 %r					- 1e <r>
#define SVMOPCODE_PUSH32_R		0x1e

// vmexit						- 1f
#define SVMOPCODE_VMEXIT		0x1f

// movw r, disp(%rip)			- 20 <r> <d> <d> <d> <d>
#define SVMOPCODE_MOVW_R_RIP	0x20
// movw r, disp8(%rbp)			- 21 <r> <d>
#define SVMOPCODE_MOVW_R_RBP8	0x21
// movw r, disp(%rbp)			- 22 <r> <d> <d> <d> <d>
#define SVMOPCODE_MOVW_R_RBP	0x22
// movw %s, (%d)				- 23 <s> <d>
#define SVMOPCODE_MOVW_R_IND	0x23
// movswl %s, %d				- 24 <s> <d>
#define SVMOPCODE_MOVSWL_R_R	0x24
// movswl disp(%rip), %r		- 25 <d> <d> <d> <d> <r>
#define SVMOPCODE_MOVSWL_RIP_R	0x25
// movswl disp(%rbp), r			- 26 <d> <d> <d> <d> <r>
#define SVMOPCODE_MOVSWL_RBP_R	0x26
// movswl (%s), %d				- 27 <s> <d>
#define SVMOPCODE_MOVSWL_IND_R	0x27
// movzwl %s, %d				- 28 <s> <d>
#define SVMOPCODE_MOVZWL_R_R	0x28
// movzwl disp(%rip), %r		- 29 <d> <d> <d> <d> <r>
#define SVMOPCODE_MOVZWL_RIP_R	0x29
// movzwl disp(%rbp), r			- 2a <d> <d> <d> <d> <r>
#define SVMOPCODE_MOVZWL_RBP_R	0x2a
// movzwl (%s), %d				- 2b <s> <d>
#define SVMOPCODE_MOVZWL_IND_R	0x2b

// ldsp %r						- 2c <r>
#define SVMOPCODE_LDSP_R		0x2c
// ldsp $disp(%rbp)				- 2d <d> <d> <d> <d>
#define SVMOPCODE_LDSP_RBP		0x2d
// ldsp (%r)					- 2e <r>
#define SVMOPCODE_LDSP_IND		0x2e

// ldbp %r						- 2f <r>
#define SVMOPCODE_LDBP_R		0x2f

// movl r, disp(%rip)			- 30 <r> <d> <d> <d> <d>
#define SVMOPCODE_MOVL_R_RIP	0x30
// movl r, disp8(%rbp)			- 31 <r> <d>
#define SVMOPCODE_MOVL_R_RBP8	0x31
// movl r, disp(%rbp)			- 32 <r> <d> <d> <d> <d>
#define SVMOPCODE_MOVL_R_RBP	0x32
// movl %s, (%d)				- 33 <s> <d>
#define SVMOPCODE_MOVL_R_IND	0x33
// movslq %s, %d				- 34 <s> <d>
#define SVMOPCODE_MOVSLQ_R_R	0x34
// movslq disp(%rip), %r		- 35 <d> <d> <d> <d> <r>
#define SVMOPCODE_MOVSLQ_RIP_R	0x35
// movslq disp(%rbp), r			- 36 <d> <d> <d> <d> <r>
#define SVMOPCODE_MOVSLQ_RBP_R	0x36
// movslq (%s), %d				- 37 <s> <d>
#define SVMOPCODE_MOVSLQ_IND_R	0x37
// movzlq %s, %d				- 38 <s> <d>
#define SVMOPCODE_MOVZLQ_R_R	0x38
// movzlq disp(%rip), %r		- 39 <d> <d> <d> <d> <r>
#define SVMOPCODE_MOVZLQ_RIP_R	0x39
// movzlq disp(%rbp), r			- 3a <d> <d> <d> <d> <r>
#define SVMOPCODE_MOVZLQ_RBP_R	0x3a
// movzlq (%s), %d				- 3b <s> <d>
#define SVMOPCODE_MOVZLQ_IND_R	0x3b

// stsp %r						- 3c <r>
#define SVMOPCODE_STSP_R		0x3c
// stsp $disp(%rbp)				- 3d <d> <d> <d> <d>
#define SVMOPCODE_STSP_RBP		0x3d
// stsp (%r)					- 3e <r>
#define SVMOPCODE_STSP_IND		0x3e

// stbp %r						- 3f <r>
#define SVMOPCODE_STBP_R		0x3f

// movq r, disp(%rip)			- 40 <r> <d> <d> <d> <d>
#define SVMOPCODE_MOVQ_R_RIP	0x40
// movq r, disp8(%rbp)			- 41 <r> <d>
#define SVMOPCODE_MOVQ_R_RBP8	0x41
// movq r, disp(%rbp)			- 42 <r> <d> <d> <d> <d>
#define SVMOPCODE_MOVQ_R_RBP	0x42
// movq %s, (%d)				- 43 <s> <d>
#define SVMOPCODE_MOVQ_R_IND	0x43

// add $imm8, %rsp				- 44 <i>
#define SVMOPCODE_ADDSP_IM8S	0x44
// add $imm, %rsp				- 45 <i> <i> <i> <i>
#define SVMOPCODE_ADDSP_IMM		0x45
// sub %r, %rsp					- 46 <r>
#define SVMOPCODE_SUBSP_R		0x46
// and $imm8, %rsp				- 47 <i>
#define SVMOPCODE_ANDSP_IM8S	0x47

// enter $stacksize				- 48 <i> <i> <i> <i>
#define SVMOPCODE_ENTER			0x48
// leave						- 49
#define SVMOPCODE_LEAVE			0x49
// ret							- 4a
#define SVMOPCODE_RET			0x4a
// retn $imm16					- 4b <i> <i>
#define SVMOPCODE_RETN			0x4b

// ncall0 %r					- 4c <r>
#define SVMOPCODE_NCALL0_R		0x4c
// ncall $args_size, %r			- 4d <a> <a> <r>
#define SVMOPCODE_NCALL_IMM_R	0x4d
// ntcall $args_size, %r		- 4e <a> <a> <r>
#define SVMOPCODE_NTCALL_IMM_R	0x4e

// nop							- 4f
#define SVMOPCODE_NOP			0x4f

// setc r						- 50 <r>
// setb r						- 50 <r>
// setnae r						- 50 <r>
#define SVMOPCODE_SETC			0x50
#define SVMOPCODE_SETB			0x50
#define SVMOPCODE_SETNAE		0x50
// setnc r						- 51 <r>
// setnb r						- 51 <r>
// setae r						- 51 <r>
#define SVMOPCODE_SETNC			0x51
#define SVMOPCODE_SETNB			0x51
#define SVMOPCODE_SETAE			0x51
// setz r						- 52 <r>
// sete r						- 52 <r>
#define SVMOPCODE_SETZ			0x52
#define SVMOPCODE_SETE			0x52
// setnz r						- 53 <r>
// setne r						- 53 <r>
#define SVMOPCODE_SETNZ			0x53
#define SVMOPCODE_SETNE			0x53
// setbe r						- 54 <r>
// setna r						- 54 <r>
#define SVMOPCODE_SETBE			0x54
#define SVMOPCODE_SETNA			0x54
// setnbe r						- 55 <r>
// seta r						- 55 <r>
#define SVMOPCODE_SETNBE		0x55
#define SVMOPCODE_SETA			0x55
// sets r						- 56 <r>
#define SVMOPCODE_SETS			0x56
// setns r						- 57 <r>
#define SVMOPCODE_SETNS			0x57
// setl r						- 58 <r>
// setnge r						- 58 <r>
#define SVMOPCODE_SETL			0x58
#define SVMOPCODE_SETNGE		0x58
// setnl r						- 59 <r>
// setge r						- 59 <r>
#define SVMOPCODE_SETNL			0x59
#define SVMOPCODE_SETGE			0x59
// setng r						- 5a <r>
// setle r						- 5a <r>
#define SVMOPCODE_SETNG			0x5a
#define SVMOPCODE_SETLE			0x5a
// setg r						- 5b <r>
// setnle r						- 5b <r>
#define SVMOPCODE_SETG			0x5b
#define SVMOPCODE_SETNLE		0x5b
// setp r						- 5c <r>
// setpe r						- 5c <r>
#define SVMOPCODE_SETP			0x5c
// setnp r						- 5d <r>
// setpo r						- 5c <r>
#define SVMOPCODE_SETNP			0x5d
// seto r						- 5e <r>
#define SVMOPCODE_SETO			0x5e
// setno r						- 5f <r>
#define SVMOPCODE_SETNO			0x5f

// jc $rel						- 60 <i> <i> <i> <i>
// jb $rel						- 60 <i> <i> <i> <i>
// jnae $rel					- 60 <i> <i> <i> <i>
#define SVMOPCODE_JC			0x60
#define SVMOPCODE_JB			0x60
#define SVMOPCODE_JNAE			0x60
// jnc $rel						- 61 <i> <i> <i> <i>
// jnb $rel						- 61 <i> <i> <i> <i>
// jae $rel						- 61 <i> <i> <i> <i>
#define SVMOPCODE_JNC			0x61
#define SVMOPCODE_JNB			0x61
#define SVMOPCODE_JAE			0x61
// jz $rel						- 62 <i> <i> <i> <i>
// je $rel						- 62 <i> <i> <i> <i>
#define SVMOPCODE_JZ			0x62
#define SVMOPCODE_JE			0x62
// jnz $rel						- 63 <i> <i> <i> <i>
// jne $rel						- 63 <i> <i> <i> <i>
#define SVMOPCODE_JNZ			0x63
#define SVMOPCODE_JNE			0x63
// jbe $rel						- 64 <i> <i> <i> <i>
// jna $rel						- 64 <i> <i> <i> <i>
#define SVMOPCODE_JBE			0x64
#define SVMOPCODE_JNA			0x64
// jnbe $rel					- 65 <i> <i> <i> <i>
// ja $rel						- 65 <i> <i> <i> <i>
#define SVMOPCODE_JNBE			0x65
#define SVMOPCODE_JA			0x65
// js $rel						- 66 <i> <i> <i> <i>
#define SVMOPCODE_JS			0x66
// jns $rel						- 67 <i> <i> <i> <i>
#define SVMOPCODE_JNS			0x67
// jl $rel						- 68 <i> <i> <i> <i>
// jnge $rel					- 68 <i> <i> <i> <i>
#define SVMOPCODE_JL			0x68
#define SVMOPCODE_JNGE			0x68
// jnl $rel						- 69 <i> <i> <i> <i>
// jge $rel						- 69 <i> <i> <i> <i>
#define SVMOPCODE_JNL			0x69
#define SVMOPCODE_JGE			0x69
// jng $rel						- 6a <i> <i> <i> <i>
// jle $rel						- 6a <i> <i> <i> <i>
#define SVMOPCODE_JNG			0x6a
#define SVMOPCODE_JLE			0x6a
// jg $rel						- 6b <i> <i> <i> <i>
// jnle $rel					- 6b  <i> <i> <i> <i>
#define SVMOPCODE_JG			0x6b
#define SVMOPCODE_JNLE			0x6b
// jo $rel						- 6c <i> <i> <i> <i>
#define SVMOPCODE_JO			0x6c
// jno $rel						- 6d <i> <i> <i> <i>
#define SVMOPCODE_JNO			0x6d
// jp $rel						- 6e <i> <i> <i> <i>
#define SVMOPCODE_JP			0x6e
// jnp $rel						- 6f <i> <i> <i> <i>
#define SVMOPCODE_JNP			0x6f

// shll $r, %r					- 70 <r> <r>
#define SVMOPCODE_SHLL_R_R		0x70
// shrl $r, %r					- 71 <r> <r>
#define SVMOPCODE_SHRL_R_R		0x71
// sarl $r, %r					- 72 <r> <r>
#define SVMOPCODE_SARL_R_R		0x72

// jmp  %r						- 73 <r>
#define SVMOPCODE_JMP_R			0x73

// shll $imm8, %r				- 74 <i> <r>
#define SVMOPCODE_SHLL_IMM_R	0x74
// shrl $imm8, %r				- 75 <i> <r>
#define SVMOPCODE_SHRL_IMM_R	0x75
// sarl $imm8, %r				- 76 <i> <r>
#define SVMOPCODE_SARL_IMM_R	0x76

// jmps $rel8					- 77 <l>
#define SVMOPCODE_JMPS			0x77

// shlq $r, %r					- 78 <r> <r>
#define SVMOPCODE_SHLQ_R_R		0x78
// shrq $r, %r					- 79 <r> <r>
#define SVMOPCODE_SHRQ_R_R		0x79
// sarq $r, %r					- 7a <r> <r>
#define SVMOPCODE_SARQ_R_R		0x7a

// jmp  $rel					- 7b <l> <l> <l> <l>
#define SVMOPCODE_JMP			0x7b

// shlq $imm8, %r				- 7c <i> <r>
#define SVMOPCODE_SHLQ_IMM_R	0x7c
// shrq $imm8, %r				- 7d <i> <r>
#define SVMOPCODE_SHRQ_IMM_R	0x7d
// sarq $imm8, %r				- 7e <i> <r>
#define SVMOPCODE_SARQ_IMM_R	0x7e

// jmp  *ind					- 7f <i> <i> <i> <i>
#define SVMOPCODE_JMP_IND		0x7f

// addl $imm8, %r				- 80 <i> <r>
#define SVMOPCODE_ADDL_IM8S_R	0x80
// addl $imm32, %r				- 81 <i> <i> <i> <i> <r>
#define SVMOPCODE_ADDL_IMM_R	0x81
// addq $imm8, %r				- 82 <i> <r>
#define SVMOPCODE_ADDQ_IM8S_R	0x82
// addq $imm64, %r				- 83 <i> <i> <i> <i> <i> <i> <i> <i> <r>
#define SVMOPCODE_ADDQ_IMM_R	0x83

// orl $imm8, %r				- 84 <i> <r>
#define SVMOPCODE_ORL_IM8S_R	0x84
// orl $imm32, %r				- 85 <i> <i> <i> <i> <r>
#define SVMOPCODE_ORL_IMM_R		0x85
// orq $imm8, %r				- 86 <i> <r>
#define SVMOPCODE_ORQ_IM8S_R	0x86
// orq $imm64, %r				- 87 <i> <i> <i> <i> <i> <i> <i> <i> <r>
#define SVMOPCODE_ORQ_IMM_R		0x87

// adcl $imm8, %r				- 88 <i> <r>
#define SVMOPCODE_ADCL_IM8S_R	0x88
// adcl $imm32, %r				- 89 <i> <i> <i> <i> <r>
#define SVMOPCODE_ADCL_IMM_R	0x89
// adcq $imm8, %r				- 8a <i> <r>
#define SVMOPCODE_ADCQ_IM8S_R	0x8a
// adcq $imm64, %r				- 8b <i> <i> <i> <i> <i> <i> <i> <i> <r>
#define SVMOPCODE_ADCQ_IMM_R	0x8b

// sbbl $imm8, %r				- 8c <i> <r>
#define SVMOPCODE_SBBL_IM8S_R	0x8c
// sbbl $imm32, %r				- 8d <i> <i> <i> <i> <r>
#define SVMOPCODE_SBBL_IMM_R	0x8d
// sbbq $imm8, %r				- 8e <i> <r>
#define SVMOPCODE_SBBQ_IM8S_R	0x8e
// sbbq $imm64, %r				- 8f <i> <i> <i> <i> <i> <i> <i> <i> <r>
#define SVMOPCODE_SBBQ_IMM_R	0x8f

// andl $imm8, %r				- 90 <i> <r>
#define SVMOPCODE_ANDL_IM8S_R	0x90
// andl $imm32, %r				- 91 <i> <i> <i> <i> <r>
#define SVMOPCODE_ANDL_IMM_R	0x91
// andq $imm8, %r				- 92 <i> <r>
#define SVMOPCODE_ANDQ_IM8S_R	0x92
// andq $imm64, %r				- 93 <i> <i> <i> <i> <i> <i> <i> <i> <r>
#define SVMOPCODE_ANDQ_IMM_R	0x93

// subl $imm8, %r				- 94 <i> <r>
#define SVMOPCODE_SUBL_IM8S_R	0x94
// subl $imm32, %r				- 95 <i> <i> <i> <i> <r>
#define SVMOPCODE_SUBL_IMM_R	0x95
// subq $imm8, %r				- 96 <i> <r>
#define SVMOPCODE_SUBQ_IM8S_R	0x96
// subq $imm64, %r				- 97 <i> <i> <i> <i> <i> <i> <i> <i> <r>
#define SVMOPCODE_SUBQ_IMM_R	0x97

// xorl $imm8, %r				- 98 <i> <r>
#define SVMOPCODE_XORL_IM8S_R	0x98
// xorl $imm32, %r				- 99 <i> <i> <i> <i> <r>
#define SVMOPCODE_XORL_IMM_R	0x99
// xorq $imm8, %r				- 9a <i> <r>
#define SVMOPCODE_XORQ_IM8S_R	0x9a
// xorq $imm64, %r				- 9b <i> <i> <i> <i> <i> <i> <i> <i> <r>
#define SVMOPCODE_XORQ_IMM_R	0x9b

// cmpl $imm8, %r				- 9c <i> <r>
#define SVMOPCODE_CMPL_IM8S_R	0x9c
// cmpl $imm32, %r				- 9d <i> <i> <i> <i> <r>
#define SVMOPCODE_CMPL_IMM_R	0x9d
// cmpq $imm8, %r				- 9e <i> <r>
#define SVMOPCODE_CMPQ_IM8S_R	0x9e
// cmpq $imm64, %r				- 9f <i> <i> <i> <i> <i> <i> <i> <i> <r>
#define SVMOPCODE_CMPQ_IMM_R	0x9f

// addl %r, %r					- a0 <r> <r>
#define SVMOPCODE_ADDL_R_R		0xa0
// addq %r, %r					- a1 <r> <r>
#define SVMOPCODE_ADDQ_R_R		0xa1

// orl %r, %r					- a2 <r> <r>
#define SVMOPCODE_ORL_R_R		0xa2
// orq %r, %r					- a3 <r> <r>
#define SVMOPCODE_ORQ_R_R		0xa3

// adcl %r, %r					- a4 <r> <r>
#define SVMOPCODE_ADCL_R_R		0xa4
// adcq %r, %r					- a5 <r> <r>
#define SVMOPCODE_ADCQ_R_R		0xa5

// sbbl %r, %r					- a6 <r> <r>
#define SVMOPCODE_SBBL_R_R		0xa6
// sbbq %r, %r					- a7 <r> <r>
#define SVMOPCODE_SBBQ_R_R		0xa7

// andl %r, %r					- a8 <r> <r>
#define SVMOPCODE_ANDL_R_R		0xa8
// andq %r, %r					- a9 <r> <r>
#define SVMOPCODE_ANDQ_R_R		0xa9

// subl %r, %r					- aa <r> <r>
#define SVMOPCODE_SUBL_R_R		0xaa
// subq %r, %r					- ab <r> <r>
#define SVMOPCODE_SUBQ_R_R		0xab

// xorl %r, %r					- ac <r> <r>
#define SVMOPCODE_XORL_R_R		0xac
// xorq %r, %r					- ad <r> <r>
#define SVMOPCODE_XORQ_R_R		0xad

// cmpl %r, %r					- ae <r> <r>
#define SVMOPCODE_CMPL_R_R		0xae
// cmpq %r, %r					- af <r> <r>
#define SVMOPCODE_CMPQ_R_R		0xaf

// imull %r1, %r2				- b0 <r1> <r2>
#define SVMOPCODE_IMULL_R_R		0xb0
// imulq %r1, %r2				- b1 <r1> <r2>
#define SVMOPCODE_IMULQ_R_R		0xb1

// divl %r1, %r2				- b2 <r1> <r2>
#define SVMOPCODE_DIVL_R_R		0xb2
// divq %r1, %r2				- b3 <r1> <r2>
#define SVMOPCODE_DIVQ_R_R		0xb3

// udivl %r1, %r2				- b4 <r1> <r2>
#define SVMOPCODE_UDIVL_R_R		0xb4
// udivq %r1, %r2				- b5 <r1> <r2>
#define SVMOPCODE_UDIVQ_R_R		0xb5

// modl %r1, %r2				- b6 <r1> <r2>
#define SVMOPCODE_MODL_R_R		0xb6
// modq %r1, %r2				- b7 <r1> <r2>
#define SVMOPCODE_MODQ_R_R		0xb7

// umodl %r1, %r2				- b8 <r1> <r2>
#define SVMOPCODE_UMODL_R_R		0xb8
// umodq %r1, %r2				- b9 <r1> <r2>
#define SVMOPCODE_UMODQ_R_R		0xb9

// fcmps %r, %r					- ba <r> <r>
#define SVMOPCODE_FCMPS_R_R		0xba
// fcmpd %r, %r					- bb <r> <r>
#define SVMOPCODE_FCMPD_R_R		0xbb
// fucmps %r, %r				- bc <r> <r>
#define SVMOPCODE_FUCMPS_R_R	0xbc
// fucmpd %r, %r				- bd <r> <r>
#define SVMOPCODE_FUCMPD_R_R	0xbd

// fadds %r, %r					- be <r> <r>
#define SVMOPCODE_FADDS_R_R		0xbe
// faddd %r, %r					- bf <r> <r>
#define SVMOPCODE_FADDD_R_R		0xbf

// fsubs %r, %r					- c0 <r> <r>
#define SVMOPCODE_FSUBS_R_R		0xc0
// fsubd %r, %r					- c1 <r> <r>
#define SVMOPCODE_FSUBD_R_R		0xc1

// fmuls %r, %r					- c2 <r> <r>
#define SVMOPCODE_FMULS_R_R		0xc2
// fmuld %r, %r					- c3 <r> <r>
#define SVMOPCODE_FMULD_R_R		0xc3

// fdivs %r, %r					- c4 <r> <r>
#define SVMOPCODE_FDIVS_R_R		0xc4
// fdivd %r, %r					- c5 <r> <r>
#define SVMOPCODE_FDIVD_R_R		0xc5

// fcvtl2s %r, %r				- c6 <r> <r>
#define SVMOPCODE_FCVTL2S_R_R	0xc6
// fcvtl2d %r, %r				- c7 <r> <r>
#define SVMOPCODE_FCVTL2D_R_R	0xc7

// fcvtq2s %r, %r				- c8 <r> <r>
#define SVMOPCODE_FCVTQ2S_R_R	0xc8
// fcvtq2d %r, %r				- c9 <r> <r>
#define SVMOPCODE_FCVTQ2D_R_R	0xc9

// fcvts2l %r, %r				- ca <r> <r>
#define SVMOPCODE_FCVTS2L_R_R	0xca
// fcvtd2l %r, %r				- cb <r> <r>
#define SVMOPCODE_FCVTD2L_R_R	0xcb

// fcvts2q %r, %r				- cc <r> <r>
#define SVMOPCODE_FCVTS2Q_R_R	0xcc
// fcvtd2q %r, %r				- cd <r> <r>
#define SVMOPCODE_FCVTD2Q_R_R	0xcd

// fcvtuq2s %r, %r				- ce <r> <r>
#define SVMOPCODE_FCVTUQ2S_R_R	0xce
// fcvtuq2d %r, %r				- cf <r> <r>
#define SVMOPCODE_FCVTUQ2D_R_R	0xcf

// fcvts2uq %r, %r				- d0 <r> <r>
#define SVMOPCODE_FCVTS2UQ_R_R	0xd0
// fcvtd2uq %r, %r				- d1 <r> <r>
#define SVMOPCODE_FCVTD2UQ_R_R	0xd1

// fcvts2d %r, %r				- d2 <r> <r>
#define SVMOPCODE_FCVTS2D_R_R	0xd2
// fcvtd2s %r, %r				- d3 <r> <r>
#define SVMOPCODE_FCVTD2S_R_R	0xd3

// lea disp(%rip), r			- d4 <d> <d> <d> <d> <r>
#define SVMOPCODE_LEA_RIP_R		0xd4
// lea disp(%rbp), r			- d5 <d> <d> <d> <d> <r>
#define SVMOPCODE_LEA_RBP_R		0xd5

// dbgbreak						- ff
#define SVMOPCODE_DBGBREAK		0xff
