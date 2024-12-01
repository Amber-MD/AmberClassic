	/* Opcodes:	*/

#define	HALT		0

#define	LOD			1
#define	STO			2
#define	LDA			3
#define	IND			4
#define	IDX			5

#define	IOR			6
#define	XOR			7
#define	ANDX		8
#define	NOTX		9

#define	LES			10
#define	LEQ			11
#define	EQU			12
#define	NEQ			13
#define	GEQ			14
#define	GTR			15

#define	NEG			16
#define	ADD			17
#define	SUB			18
#define	MUL			19
#define	DIV			20
#define	MOD			21

#define	JMP			22
#define	FJP			23

#define	MRK			24
#define	CLR			25
#define	CUF			26
#define	ENT			27
#define	RET			28
#define	CGF			29
#define	CSF			30

#define	NOOP		31

#define	N_OPCODE	32


	/* System Calls:	*/

#define	SC_ADDBOND	1
#define	SC_ADDRESIDUE	2
#define	SC_ADDSTRAND	3
#define	SC_CAP		4
#define	SC_FCLOSE	5
#define	SC_FOPEN	6
#define	SC_FPRINTF	7
#define	SC_FPUTLINKIN	8
#define	SC_FPUTPDB	9
#define	SC_FSCANF 	10
#define	SC_GETRESIDUE	11
#define	SC_GETTRANSFORM	12
#define	SC_LENGTH	13
#define	SC_NEWMOLECULE  14
#define	SC_NEWTRANSFORM 15
#define	SC_PRINTF	16
#define	SC_SCANF 	17
#define	SC_SETTRANSFORM	18
#define	SC_SUBSTR	19
#define	SC_TRANSFORMRES	20
#define	SC_UPDTRANSFORM	21


	/* Machine Defs:	*/
	/* segments:		*/

#define	S_UNDEF		(-1)
#define	S_GLOBAL	0
#define	S_LOCAL		1
#define	S_LIT		2
#define	S_CODE		3
#define	S_SYSCALL	4


	/* Activation Record:	*/

#define	FR_VAL		0	/* function return val	*/
#define	R_ADDR		1	/* return address	*/
#define	O_FRMP		2	/* old frame pointer	*/
#define	O_MRKP		3	/* old mark pointer 	*/
#define	O_STKP		4	/* old stack pointer 	*/
#define	AR_SIZE		5


	/* Instruction format:	*/

typedef	struct	inst	{
	int	i_op;
	int	i_seg;
	int	i_off;
} INST;


	/* Sizes:		*/

	/* Instructions		*/

#define	IMEM_SIZE		10000

	/* Data:			*/

#define	LIT_SIZE		200
#define	LOCAL_SIZE		200
#define	GLOBAL_SIZE		200
#define	DMEM_SIZE		( GLOBAL_SIZE + LOCAL_SIZE + 1000 )
#define	DMEM_LADDR		( DMEM_SIZE - 1 )
