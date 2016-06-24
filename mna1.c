/*
Elementos aceitos e linhas do netlist:

Resistor:  R<nome> <no+> <no-> <resistencia>
VCCS:      G<nome> <io+> <io-> <vi+> <vi-> <transcondutancia>   = Fonte de tensao controlada a corrente
VCVC:      E<nome> <vo+> <vo-> <vi+> <vi-> <ganho de tensao>    = Fonte de tensao controlada a tensao
CCCS:      F<nome> <io+> <io-> <ii+> <ii-> <ganho de corrente>  = Fonte de corrente controlada a corrente
CCVS:      H<nome> <vo+> <vo-> <ii+> <ii-> <transresistencia>   = Fonte de corrente controlada a tensao
Fonte I:   I<nome> <io+> <io-> <corrente>
Fonte V:   V<nome> <vo+> <vo-> <tensao>
Amp. op.:  O<nome> <vo1> <vo2> <vi1> <vi2>

Elementos incluídos para o trabalho:

Indutor: L<nome> <nó1> <nó2> <Indutância>
Acoplamento entre indutores: K<nome> <La> <Lb> <k> (La e Lb nomes de indutores já declarados.)
Capacitor: C<nome> <nó1> <nó2> <Capacitância>
Fonte de corrente: I<nome> <nó+> <nó-> <módulo> <fase (graus)> <valor contínuo>
Fonte de tensão: V<nome> <nó+> <nó-> <módulo> <fase (graus)> <valor contínuo>
 *Transistor MOS: M<nome> <nód> <nóg> <nós> <nób> <NMOS ou PMOS> L=<comprimento> W=<largura> <K> <Vt0> <lamba> <gama> <phi> <Ld>

As fontes F e H tem o ramo de entrada em curto
O amplificador operacional ideal tem a saida suspensa
 */

#include <stdio.h>
#include <conio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include <math.h>
#include <complex.h>
#define MAX_LINHA 80
#define MAX_NOME 11
#define MAX_ELEM 50
#define MAX_NOS 50
#define TOLG 1e-9
//#define DEBUG

#define MAX_ITERACOES 50
#define MAX_INIT 5
#define ERRO_MAX 0.01
#define X_ERRO 1

typedef enum TipoPontoOperacao
{
	corte, triodo, saturacao
} TipoPontoOperacao;

typedef struct elemento // Elemento do netlist
{
	char nome[MAX_NOME], nome1[MAX_NOME], nome2[MAX_NOME];
	double valor, modulo, fase;
	double L, W, K, Vt0, lambda, gama, phi, Ld, gm, gds, gmb, i0;
	int a,b,c,d,x,y,g,s;
	TipoPontoOperacao pontoOperacao;
} elemento;


elemento netlist[MAX_ELEM]; /* Netlist */


int
numeroElementos, /* Numero de elementos */
numeroVariaveis, /* Numero de variaveis */
numeroNos, /* Numero de nos */
i,j,k,x,y,
ponto,
numeroPontos,
sistemaLinear;


char
/* Foram colocados limites nos formatos de leitura para alguma protecao
   contra excesso de caracteres nestas variaveis */
nomeArquivo[MAX_LINHA+1],
*nomeArquivoTab,
tipo,
na[MAX_NOME],nb[MAX_NOME],nc[MAX_NOME],nd[MAX_NOME],
lista[MAX_NOS+1][MAX_NOME+2], /*Tem que caber jx antes do nome (motivo do +2), imagino*/
txt[MAX_LINHA+1],
*p,
tipoAC[4]; //3+1

FILE *arquivo, *arquivoTab;

double
g,
Yn[MAX_NOS+1][MAX_NOS+2],
YnAnterior[MAX_NOS+1][MAX_NOS+2],
erros[MAX_NOS+1],
freqInicial,
freqFinal,
deltaOmega,
omega,
A, B;

double complex
gComplex,
YnComplex[MAX_NOS+1][MAX_NOS+2],
YnAnteriorComplex[MAX_NOS+1][MAX_NOS+2];

/* Resolucao de sistema de equacoes lineares.
   Metodo de Gauss-Jordan com condensacao pivotal */

elemento ProcurarIndutor (char* nome)
{
	int i;

	for (i=1; i<numeroElementos; i++) //netlist[0] nao e usado
	{
		if (strcmp (nome, netlist[i].nome) == 0)
			return netlist[i];
	}

	printf ("Indutor %s do acoplamento nao declarado\n", nome);
	exit (1);
}

int MontarCapacitanciasMOSPontoOp(elemento mos)
{
	g = 1/1e9;

	//Capacitancia G-S
	Yn[mos.s][mos.s] += g;
	Yn[mos.g][mos.g] += g;
	Yn[mos.s][mos.g] -= g;
	Yn[mos.g][mos.s] -= g;

	//Capacitancia G-D
	Yn[mos.d][mos.d] += g;
	Yn[mos.g][mos.g] += g;
	Yn[mos.d][mos.g] -= g;
	Yn[mos.g][mos.d] -= g;

	//Capacitancia G-B
	Yn[mos.b][mos.b] += g;
	Yn[mos.g][mos.g] += g;
	Yn[mos.b][mos.g] -= g;
	Yn[mos.g][mos.b] -= g;

	return 0;
}

int MontarCapacitanciasMOSAnaliseFreq(elemento mos, double w)
{
	double Cox, mi;
	double Cgs, Cgd, Cgb;

	if (strcmp(mos.nome1, "NMOS") == 0)
		mi = 0.05;
	else
		mi = 0.02;

	Cox = 2*mos.K/mi;

	if (mos.pontoOperacao == corte)
	{
		Cgs = Cox*mos.W*mos.Ld;
		Cgd = Cgs;
		Cgb = Cox*mos.W*mos.L;
	}
	else if (mos.pontoOperacao == triodo)
	{
		Cgs = Cox*mos.W*mos.Ld + Cox*mos.W*mos.L/2.0;
		Cgd = Cgs;
		Cgb = 0;
	}
	else
	{
		Cgs = Cox*mos.W*mos.Ld + 2*Cox*mos.W*mos.L/3.0;
		Cgd = Cox*mos.W*mos.Ld;
		Cgb = 0;
	}


	//Capacitancia G-S
	g = I*w*Cgs;
	Yn[mos.s][mos.s] += g;
	Yn[mos.g][mos.g] += g;
	Yn[mos.s][mos.g] -= g;
	Yn[mos.g][mos.s] -= g;

	//Capacitancia G-D
	g = I*w*Cgd;
	Yn[mos.d][mos.d] += g;
	Yn[mos.g][mos.g] += g;
	Yn[mos.d][mos.g] -= g;
	Yn[mos.g][mos.d] -= g;

	//Capacitancia G-B
	g = I*w*Cgb;
	Yn[mos.b][mos.b] += g;
	Yn[mos.g][mos.g] += g;
	Yn[mos.b][mos.g] -= g;
	Yn[mos.g][mos.b] -= g;

	return 0;
}

int MontarEstampaMOS (elemento mos)
{
	double vgs,vds,vbs,Vt;
	int auxiliar;

	if (mos.nome[0] != 'M')
		return 1;


	// Garante que D é o nó de maior tensão no NMOS
	if (strcmp (mos.nome1, "NMOS") == 0)
	{
		if (YnAnterior[mos.d][numeroVariaveis+1] < YnAnterior[mos.s][numeroVariaveis+1])
		{
			auxiliar = mos.d;
			mos.d = mos.s;
			mos.s = auxiliar;
		}
	}
	// Garante que S é o nó de maior tensão no PMOS
	else
	{
		if (YnAnterior[mos.d][numeroVariaveis+1] > YnAnterior[mos.s][numeroVariaveis+1])
		{
			auxiliar = mos.d;
			mos.d = mos.s;
			mos.s = auxiliar;
		}
	}


	vgs = YnAnterior[mos.g][numeroVariaveis+1] - YnAnterior[mos.s][numeroVariaveis+1];
	vds = YnAnterior[mos.d][numeroVariaveis+1] - YnAnterior[mos.s][numeroVariaveis+1];
	vbs = YnAnterior[mos.b][numeroVariaveis+1] - YnAnterior[mos.s][numeroVariaveis+1];

	if (vbs > mos.phi)
		vbs = mos.phi/2.0;

	Vt = mos.Vt0 + mos.gama*(sqrt(mos.phi - vbs) - sqrt(mos.phi));

	if (strcmp (netlist[numeroElementos].nome1, "PMOS") == 0)
		Vt *= -1;

	if (fabs(vgs) < fabs(Vt)) //Corte
	{
		mos.pontoOperacao = corte;
		mos.i0 = 0;
		mos.gm = 0;
		mos.gds = 0;
		mos.gmb = 0;
	}

	else if (fabs(vds) <= fabs(vgs - Vt)) //Triodo
	{
		mos.pontoOperacao = triodo;

		mos.gm = mos.K * (mos.W/mos.L) * (2*vds) * (1 + mos.lambda*vds);
		mos.gds = mos.K * (mos.W/mos.L) * (2*(vgs-Vt) - 2*vds + 4*mos.lambda*(vgs-Vt)*vds - 3*mos.lambda*pow(vds,2));
		mos.gmb = (mos.gm*mos.gama)/(2*sqrt(mos.phi-vbs));

		mos.i0 = mos.K * (mos.W/mos.L) * (2*(vgs-Vt)*vds - pow(vds,2)) * (1+mos.lambda*vds) - mos.gm*vgs - mos.gds*vds - mos.gmb*vbs;
	}

	else //Saturacao
	{
		mos.pontoOperacao = saturacao;

		//printf ("K = %f | W = %f | L = %f | Vgs = %f | Vt = %f | Lambda = %f | Vds = %f\n", mos.K,mos.W,mos.L,vgs,Vt,mos.lambda,vds);
		mos.gm = mos.K * (mos.W/mos.L) * 2 * (vgs-Vt) * (1+mos.lambda*vds);
		mos.gds = mos.K * (mos.W/mos.L) * mos.lambda * pow (vgs-Vt,2);
		mos.gmb = (mos.gm*mos.gama)/(2*sqrt(mos.phi-vbs));

		mos.i0 = mos.K * (mos.W/mos.L) * pow(vgs-Vt, 2) * (1+mos.lambda*vds) - mos.gm*vgs - mos.gds*vds - mos.gmb*vbs;
	}

	Yn[mos.d][mos.d] += mos.gds;
	Yn[mos.d][mos.s] -= mos.gds + mos.gm + mos.gmb;
    Yn[mos.d][mos.g] += mos.gm;
    Yn[mos.d][mos.b] += mos.gmb;
    Yn[mos.s][mos.d] -= mos.gds;
    Yn[mos.s][mos.s] += mos.gds + mos.gm + mos.gmb;
    Yn[mos.s][mos.g] -= mos.gm;
    Yn[mos.s][mos.b] -= mos.gmb;

    Yn[mos.d][numeroVariaveis+1] -= mos.i0;
    Yn[mos.s][numeroVariaveis+1] += mos.i0;

    printf ("\n");
    switch (mos.pontoOperacao)
    {
    	case corte:
    		printf ("Modo = Corte\n");
    		break;
    	case triodo:
    		printf ("Modo = Triodo\n");
    		break;
    	case saturacao:
    		printf ("Modo = Saturacao\n");
    		break;
    }

    printf ("\n");
    //printf ("d = %i | g = %i | s = %i | b = %i\n", mos.d, mos.g, mos.s, mos.b);

    printf ("Vd = %f | Vg = %f | Vs = %f | Vb = %f\n", YnAnterior[mos.d][numeroVariaveis+1],
    												   YnAnterior[mos.g][numeroVariaveis+1],
													   YnAnterior[mos.s][numeroVariaveis+1],
													   YnAnterior[mos.b][numeroVariaveis+1]);

    printf ("I0 = %f | Gm = %f | Gds = %f | Gmb = %f | Vt = %f\n", mos.i0, mos.gm, mos.gds, mos.gmb, Vt);

	return 0;
}

int ResolverSistema(void)
{
	int i,j,l, a;
	double t, p;

	for (i=1; i<=numeroVariaveis; i++) {
		t=0.0;
		a=i;
		for (l=i; l<=numeroVariaveis; l++) {
			if (fabs(Yn[l][i])>fabs(t)) {
				a=l;
				t=Yn[l][i];
			}
		}
		if (i!=a) {
			for (l=1; l<=numeroVariaveis+1; l++) {
				p=Yn[i][l];
				Yn[i][l]=Yn[a][l];
				Yn[a][l]=p;
			}
		}
		if (fabs(t)<TOLG) {
			printf("Sistema singular\n");
			return 1;
		}
		for (j=numeroVariaveis+1; j>0; j--) {  /* Basta j>i em vez de j>0 */
			Yn[i][j]/= t;
			p=Yn[i][j];
			if (p!=0)  /* Evita operacoes com zero */
				for (l=1; l<=numeroVariaveis; l++) {
					if (l!=i)
						Yn[l][j]-=Yn[l][i]*p;
				}
		}
	}
	return 0;
}

int ResolverSistemaComplexo(void)
{
	int i,j,l, a;
	double complex t, p;

	for (i=1; i<=numeroVariaveis; i++) {
		t=0.0;
		a=i;
		for (l=i; l<=numeroVariaveis; l++) {
			if (cabs(YnComplex[l][i])>cabs(t)) {
				a=l;
				t=YnComplex[l][i];
			}
		}
		if (i!=a) {
			for (l=1; l<=numeroVariaveis+1; l++) {
				p=YnComplex[i][l];
				YnComplex[i][l]=YnComplex[a][l];
				YnComplex[a][l]=p;
			}
		}
		if (cabs(t)<TOLG) {
			printf("Sistema singular\n");
			return 1;
		}
		for (j=numeroVariaveis+1; j>0; j--) {  /* Basta j>i em vez de j>0 */
			YnComplex[i][j]/= t;
			p=YnComplex[i][j];
			if (p!=0)  /* Evita operacoes com zero */
				for (l=1; l<=numeroVariaveis; l++)
				{
					if (l!=i)
						YnComplex[l][j]-=YnComplex[l][i]*p;
				}
		}
	}
	return 0;
}

void SalvarUltimaIteracao()
{
	for (i=1; i<=numeroVariaveis; i++)
	{
		for (j=1; j<=numeroVariaveis+1; j++)
			YnAnterior[i][j] = Yn[i][j];
	}
}

/* Rotina que conta os nos e atribui numeros a eles */
int NumerarNo(char *nome)
{
	int i,achou;

	i=0; achou=0;
	while (!achou && i<=numeroVariaveis)  // Procura pra ver se o nó "nome" já existe. Enquanto não acha, vai contando quantos existem
		if (!(achou=!strcmp(nome,lista[i])))
			i++;

	if (!achou)
	{
		if (numeroVariaveis==MAX_NOS)
		{
			printf("O programa so aceita ate %d nos\n",numeroVariaveis);
			exit(1);
		}
		numeroVariaveis++;
		strcpy(lista[numeroVariaveis],nome);     // Adiciona o novo nó ao array de nós
		return numeroVariaveis; /* novo no */
	}
	else
	{
		return i; /* no ja conhecido */
	}
}

void ZerarSistema()
{
	for (i=0; i<=numeroVariaveis; i++)
	{
		for (j=0; j<=numeroVariaveis+1; j++)
			Yn[i][j]=0;
	}
}

void ZerarSistemaComplexo()
{
	for (i=0; i<=numeroVariaveis; i++)
	{
		for (j=0; j<=numeroVariaveis+1; j++)
			YnComplex[i][j]=0;
	}
}

void PrintarSistema(double sistema[MAX_NOS+1][MAX_NOS+2])
{
	for (i=1; i<=numeroVariaveis; i++)
	{
		for (j=1; j<=numeroVariaveis+1; j++)
			if (sistema[i][j]!=0) printf("%f ",sistema[i][j]);
			else printf(" ... ");
		printf("\n");
	}
}

void PrintarSistemaComplexo()
{
	for (k=1; k<=numeroVariaveis; k++)
	{
		for (j=1; j<=numeroVariaveis+1; j++)
			if (cabs(YnComplex[k][j]) != 0)
				printf("%3.1f+%3.1fj ", creal(YnComplex[k][j]), cimag(YnComplex[k][j]));
			else printf(" ... ");
		printf("\n");
	}
}

void MontarEstampasPontoOp()
{
	for (i=1; i<=numeroElementos; i++)
	{
		tipo=netlist[i].nome[0];
		if (tipo=='R')
		{
			g = 1/netlist[i].valor;
			Yn[netlist[i].a][netlist[i].a] += g;
			Yn[netlist[i].b][netlist[i].b] += g;
			Yn[netlist[i].a][netlist[i].b] -= g;
			Yn[netlist[i].b][netlist[i].a] -= g;
		}
		else if (tipo=='L')
		{
			// Vira um R de 1nOhms
			g = 1/1e-9;
			Yn[netlist[i].a][netlist[i].a] += g;
			Yn[netlist[i].b][netlist[i].b] += g;
			Yn[netlist[i].a][netlist[i].b] -= g;
			Yn[netlist[i].b][netlist[i].a] -= g;
		}
		else if (tipo=='C')
		{
			// Vira um R de 1GOhms
			g = 1/1e9;
			Yn[netlist[i].a][netlist[i].a] += g;
			Yn[netlist[i].b][netlist[i].b] += g;
			Yn[netlist[i].a][netlist[i].b] -= g;
			Yn[netlist[i].b][netlist[i].a] -= g;
		}
		else if (tipo=='G')
		{
			g = netlist[i].valor;
			Yn[netlist[i].a][netlist[i].c] += g;
			Yn[netlist[i].b][netlist[i].d] += g;
			Yn[netlist[i].a][netlist[i].d] -= g;
			Yn[netlist[i].b][netlist[i].c] -= g;
		}
		else if (tipo=='I')  // I e V só o valor contínuo
		{
			g = netlist[i].valor;
			Yn[netlist[i].a][numeroVariaveis+1] -= g;
			Yn[netlist[i].b][numeroVariaveis+1] += g;
		}
		else if (tipo=='V')
		{
			Yn[netlist[i].a][netlist[i].x] += 1;
			Yn[netlist[i].b][netlist[i].x] -= 1;
			Yn[netlist[i].x][netlist[i].a] -= 1;
			Yn[netlist[i].x][netlist[i].b] += 1;
			Yn[netlist[i].x][numeroVariaveis+1] -= netlist[i].valor;
		}
		else if (tipo=='E')
		{
			g = netlist[i].valor;
			Yn[netlist[i].a][netlist[i].x] += 1;
			Yn[netlist[i].b][netlist[i].x] -= 1;
			Yn[netlist[i].x][netlist[i].a] -= 1;
			Yn[netlist[i].x][netlist[i].b] += 1;
			Yn[netlist[i].x][netlist[i].c] += g;
			Yn[netlist[i].x][netlist[i].d] -= g;
		}
		else if (tipo=='F')
		{
			g = netlist[i].valor;
			Yn[netlist[i].a][netlist[i].x] += g;
			Yn[netlist[i].b][netlist[i].x] -= g;
			Yn[netlist[i].c][netlist[i].x] += 1;
			Yn[netlist[i].d][netlist[i].x] -= 1;
			Yn[netlist[i].x][netlist[i].c] -= 1;
			Yn[netlist[i].x][netlist[i].d] += 1;
		}
		else if (tipo=='H')
		{
			g = netlist[i].valor;
			Yn[netlist[i].a][netlist[i].y] += 1;
			Yn[netlist[i].b][netlist[i].y] -= 1;
			Yn[netlist[i].c][netlist[i].x] += 1;
			Yn[netlist[i].d][netlist[i].x] -= 1;
			Yn[netlist[i].y][netlist[i].a] -= 1;
			Yn[netlist[i].y][netlist[i].b] += 1;
			Yn[netlist[i].x][netlist[i].c] -= 1;
			Yn[netlist[i].x][netlist[i].d] += 1;
			Yn[netlist[i].y][netlist[i].x] += g;
		}
		else if (tipo=='O')
		{
			Yn[netlist[i].a][netlist[i].x] += 1;
			Yn[netlist[i].b][netlist[i].x] -= 1;
			Yn[netlist[i].x][netlist[i].c] += 1;
			Yn[netlist[i].x][netlist[i].d] -= 1;
		}
		else if (tipo=='M')
		{
			MontarEstampaMOS(netlist[i]);
			MontarCapacitanciasMOSPontoOp(netlist[i]);
		}
#ifdef DEBUG
		/* Opcional: Mostra o sistema apos a montagem da estampa */
		printf("Sistema apos a estampa de %s\n",netlist[i].nome);
		for (k=1; k<=numeroVariaveis; k++)
		{
			for (j=1; j<=numeroVariaveis+1; j++)
				if (Yn[k][j]!=0)
					printf("%+3.1f ",Yn[k][j]);
				else printf(" ... ");
			printf("\n");
		}
		getch();
#endif
	}
}

// Ainda não coloquei o MOS
void MontarEstampasAnaliseFrequencia(double w)
{
	for (i=1; i<=numeroElementos; i++)
	{
		tipo=netlist[i].nome[0];
		if (tipo=='R')
		{
			gComplex = 1/netlist[i].valor;
			YnComplex[netlist[i].a][netlist[i].a] += gComplex;
			YnComplex[netlist[i].b][netlist[i].b] += gComplex;
			YnComplex[netlist[i].a][netlist[i].b] -= gComplex;
			YnComplex[netlist[i].b][netlist[i].a] -= gComplex;
		}
		else if (tipo=='L')
		{
			gComplex = 1/(I*w*netlist[i].valor); // 1/jwL
			//printf ("Modulo gL = %g\n", cabs (gComplex));
			YnComplex[netlist[i].a][netlist[i].x] += 1;
			YnComplex[netlist[i].b][netlist[i].x] -= 1;
			YnComplex[netlist[i].x][netlist[i].a] -= 1;
			YnComplex[netlist[i].x][netlist[i].b] += 1;
			YnComplex[netlist[i].x][netlist[i].x] += 1/gComplex;
		}
		else if (tipo=='K')
		{
			// Ainda nao testei
			// Acrescenta a parte relativa ao acoplamento j*w*M
			gComplex = I*w*netlist[i].valor;
			x = ProcurarIndutor (netlist[i].nome1).x;
			y = ProcurarIndutor (netlist[i].nome2).x;
			YnComplex[x][y] += gComplex;
			YnComplex[y][x] += gComplex;
		}
		else if (tipo=='C')
		{
			gComplex = I*w*netlist[i].valor;  // jwC
			//printf ("Capacitor: %g + %gi\n", creal(gComplex), cimag(gComplex));
			YnComplex[netlist[i].a][netlist[i].a] += gComplex;
			YnComplex[netlist[i].b][netlist[i].b] += gComplex;
			YnComplex[netlist[i].a][netlist[i].b] -= gComplex;
			YnComplex[netlist[i].b][netlist[i].a] -= gComplex;
		}
		else if (tipo=='G')
		{
			gComplex = netlist[i].valor;
			YnComplex[netlist[i].a][netlist[i].c] += gComplex;
			YnComplex[netlist[i].b][netlist[i].d] += gComplex;
			YnComplex[netlist[i].a][netlist[i].d] -= gComplex;
			YnComplex[netlist[i].b][netlist[i].c] -= gComplex;
		}
		else if (tipo=='I')
		{
			A = cos(netlist[i].fase*M_PI/180.0) * netlist[i].modulo;
			B = sin(netlist[i].fase*M_PI/180.0) * netlist[i].modulo;
			gComplex = A+B*I;

			YnComplex[netlist[i].a][numeroVariaveis+1] -= gComplex;
			YnComplex[netlist[i].b][numeroVariaveis+1] += gComplex;
		}
		else if (tipo=='V')
		{
			A = cos(netlist[i].fase*M_PI/180.0) * netlist[i].modulo;
			B = sin(netlist[i].fase*M_PI/180.0) * netlist[i].modulo;
			gComplex = A+B*I;

			YnComplex[netlist[i].a][netlist[i].x] += 1;
			YnComplex[netlist[i].b][netlist[i].x] -= 1;
			YnComplex[netlist[i].x][netlist[i].a] -= 1;
			YnComplex[netlist[i].x][netlist[i].b] += 1;
			YnComplex[netlist[i].x][numeroVariaveis+1] -= gComplex;
		}
		else if (tipo=='E')
		{
			gComplex = netlist[i].valor;
			YnComplex[netlist[i].a][netlist[i].x] += 1;
			YnComplex[netlist[i].b][netlist[i].x] -= 1;
			YnComplex[netlist[i].x][netlist[i].a] -= 1;
			YnComplex[netlist[i].x][netlist[i].b] += 1;
			YnComplex[netlist[i].x][netlist[i].c] += gComplex;
			YnComplex[netlist[i].x][netlist[i].d] -= gComplex;
		}
		else if (tipo=='F')
		{
			gComplex = netlist[i].valor;
			YnComplex[netlist[i].a][netlist[i].x] += gComplex;
			YnComplex[netlist[i].b][netlist[i].x] -= gComplex;
			YnComplex[netlist[i].c][netlist[i].x] += 1;
			YnComplex[netlist[i].d][netlist[i].x] -= 1;
			YnComplex[netlist[i].x][netlist[i].c] -= 1;
			YnComplex[netlist[i].x][netlist[i].d] += 1;
		}
		else if (tipo=='H')
		{
			gComplex = netlist[i].valor;
			YnComplex[netlist[i].a][netlist[i].y] += 1;
			YnComplex[netlist[i].b][netlist[i].y] -= 1;
			YnComplex[netlist[i].c][netlist[i].x] += 1;
			YnComplex[netlist[i].d][netlist[i].x] -= 1;
			YnComplex[netlist[i].y][netlist[i].a] -= 1;
			YnComplex[netlist[i].y][netlist[i].b] += 1;
			YnComplex[netlist[i].x][netlist[i].c] -= 1;
			YnComplex[netlist[i].x][netlist[i].d] += 1;
			YnComplex[netlist[i].y][netlist[i].x] += gComplex;
		}
		else if (tipo=='O')
		{
			YnComplex[netlist[i].a][netlist[i].x] += 1;
			YnComplex[netlist[i].b][netlist[i].x] -= 1;
			YnComplex[netlist[i].x][netlist[i].c] += 1;
			YnComplex[netlist[i].x][netlist[i].d] -= 1;
		}
#ifdef DEBUG
		/* Opcional: Mostra o sistema apos a montagem da estampa */
		printf("Sistema apos a estampa de %s\n",netlist[i].nome);
		PrintarSistemaComplexo();
		//getch();
#endif
	}
}

int TestarConvergencia()
{
	int i;

	for (i=1; i<=numeroVariaveis;i++)
	{
		if (fabs(Yn[i][numeroVariaveis+1]) > X_ERRO)
		{
			erros[i] = X_ERRO*fabs((Yn[i][numeroVariaveis+1]-YnAnterior[i][numeroVariaveis+1])/Yn[i][numeroVariaveis+1]);
		}
		else
		{
			erros[i] = fabs(Yn[i][numeroVariaveis+1]-YnAnterior[i][numeroVariaveis+1]);
		}

		if (erros[i] > ERRO_MAX)
			return 0;
	}

	return 1;
}

void InicializacaoRandomica()
{
	double valor;

    srand( (unsigned)time(NULL) );

	for (i=1; i<=numeroVariaveis;i++)
	{
		if (erros[i] > ERRO_MAX)
		{
			valor = rand() % 201;
			valor -= 100;
			valor /= 100.0;
			YnAnterior[i][numeroVariaveis+1] = valor;
		}
		else
			YnAnterior[i][numeroVariaveis+1] = Yn[i][numeroVariaveis+1];
	}
}

void NewtonRaphsonPontoOp ()
{
	int iteracao, i, convergiu, init;

	// Zerando Yn e YnAnterior
	ZerarSistema();
	SalvarUltimaIteracao(); // Zerar o YnAnterior tambem

	//Inicializacao das variaveis do sistema
	for (i=1; i<=numeroVariaveis; i++)
	{
		YnAnterior[i][numeroVariaveis+1] = 0.1;
	}

	//Debug
	printf ("Printando YnAnterior inicializado\n");
	PrintarSistema (YnAnterior);


//	MontarEstampasPontoOp();
//	return;

	for (init=1; init<=MAX_INIT; init++)
	{
		if (init != 1)
		{
			InicializacaoRandomica();
		}

		for (iteracao = 1; iteracao<=MAX_ITERACOES; iteracao++)
		{
			ZerarSistema();
			MontarEstampasPontoOp();
			ResolverSistema();
//			if (ResolverSistema())
//			{
//				getch();
//				exit(1);
//			}

			convergiu = TestarConvergencia();
			if (convergiu == 1)
			{
				printf ("Sistema convergiu com %i iteracoes\n",iteracao);
				return;
			}
			PrintarSistema(Yn);
			SalvarUltimaIteracao();
		}
	}

	printf ("\n");
	for (i=1; i<=numeroVariaveis;i++)
	{
		printf ("Erro %i : %f\n", i, erros[i]);
	}
	printf ("Nem fodendo que essa porra converge\n");
}



int main(void)
{
	//clrscr();

	printf("Programa do trabalho de CE II 2016.1\n");
	printf("Análise de ponto de operação e de resposta em frequência de circuitos lineares contendo transistores MOS\n");
	printf("Por:\n");
	printf("   Lucas de Andrade Cerqueira\n");
	printf("   Bruno Granato\n");
	printf("   Joao Felipe Guedes\n");
	denovo:
	/* Leitura do netlist */
	sistemaLinear = 1;
	numeroElementos=0;
	numeroVariaveis=0;
	strcpy(lista[0],"0");

	printf("Nome do arquivo com o netlist (ex: mna.net): ");
	scanf("%50s",nomeArquivo);
	arquivo=fopen(nomeArquivo,"r");
	if (arquivo==0) {
		printf("Arquivo %s inexistente\n",nomeArquivo);
		goto denovo;
	}
	printf("Lendo netlist:\n");
	fgets(txt,MAX_LINHA,arquivo);
	printf("Titulo: %s",txt);
	while (fgets(txt,MAX_LINHA,arquivo))
	{
		numeroElementos++; /* Nao usa o netlist[0] */
		if (numeroElementos>MAX_ELEM)
		{
			printf("O programa so aceita ate %d elementos\n",MAX_ELEM);
			exit(1);
		}
		txt[0]=toupper(txt[0]);
		tipo=txt[0];
		sscanf(txt,"%10s",netlist[numeroElementos].nome);
		p=txt+strlen(netlist[numeroElementos].nome); /* Inicio dos parametros */
		/* O que e lido depende do tipo */
		if (tipo=='R' || tipo=='L' || tipo=='C')
		{
			sscanf(p,"%10s%10s%lg",na,nb,&netlist[numeroElementos].valor);
			printf("%s %s %s %g\n",netlist[numeroElementos].nome,na,nb,netlist[numeroElementos].valor);
			netlist[numeroElementos].a=NumerarNo(na);   //Associa na e nb do netlist aos numeros dos nos
			netlist[numeroElementos].b=NumerarNo(nb);
		}
		else if (tipo=='I' || tipo=='V')
		{
			sscanf(p,"%10s%10s%lg%lg%lg",na,nb,&netlist[numeroElementos].modulo,&netlist[numeroElementos].fase,&netlist[numeroElementos].valor);
			printf("%s %s %s %g %g %g\n",netlist[numeroElementos].nome,na,nb,
										 netlist[numeroElementos].modulo, netlist[numeroElementos].fase, netlist[numeroElementos].valor);
			netlist[numeroElementos].a=NumerarNo(na);
			netlist[numeroElementos].b=NumerarNo(nb);
		}
		else if (tipo=='G' || tipo=='E' || tipo=='F' || tipo=='H')
		{
			sscanf(p,"%10s%10s%10s%10s%lg",na,nb,nc,nd,&netlist[numeroElementos].valor);
			printf("%s %s %s %s %s %g\n",netlist[numeroElementos].nome,na,nb,nc,nd,netlist[numeroElementos].valor);
			netlist[numeroElementos].a=NumerarNo(na);
			netlist[numeroElementos].b=NumerarNo(nb);
			netlist[numeroElementos].c=NumerarNo(nc);
			netlist[numeroElementos].d=NumerarNo(nd);
		}
		else if (tipo=='O')
		{
			sscanf(p,"%10s%10s%10s%10s",na,nb,nc,nd);
			printf("%s %s %s %s %s\n",netlist[numeroElementos].nome,na,nb,nc,nd);
			netlist[numeroElementos].a=NumerarNo(na);
			netlist[numeroElementos].b=NumerarNo(nb);
			netlist[numeroElementos].c=NumerarNo(nc);
			netlist[numeroElementos].d=NumerarNo(nd);
		}

		else if (tipo=='K')
		{
			//nome1 e nome2 guardam os nomes dos indutores previamente declarados
			sscanf(p,"%10s%10s%lg",netlist[numeroElementos].nome1,netlist[numeroElementos].nome2,&netlist[numeroElementos].valor);

			double L1 = ProcurarIndutor (netlist[numeroElementos].nome1).valor;
			double L2 = ProcurarIndutor (netlist[numeroElementos].nome2).valor;

			//M = k*sqrt(L1*L2)
			netlist[numeroElementos].valor = netlist[numeroElementos].valor * sqrt (L1*L2);

			printf("%s %s %s %g\n",netlist[numeroElementos].nome,
								   netlist[numeroElementos].nome1,
								   netlist[numeroElementos].nome2,
								   netlist[numeroElementos].valor);
		}

		else if (tipo=='M')
		{
			sistemaLinear = 0;

			char L[10], W[10];
			//*Transistor MOS: M<nome> <nód> <nóg> <nós> <nób> <NMOS ou PMOS> L=<comprimento> W=<largura> <K> <Vt0> <lambda> <gama> <phi> <Ld>
			sscanf(p,"%10s%10s%10s%10s%10s%10s%10s%lg%lg%lg%lg%lg%lg",na,nb,nc,nd,netlist[numeroElementos].nome1,
										 	 	 	L, W, &netlist[numeroElementos].K,
										            &netlist[numeroElementos].Vt0, &netlist[numeroElementos].lambda, &netlist[numeroElementos].gama,
										            &netlist[numeroElementos].phi, &netlist[numeroElementos].Ld);

			char *ptr;
			char *token;
			token = strtok(L,"=");
			token= strtok(NULL,"=");
			netlist[numeroElementos].L = strtod(token, &ptr);

			token = strtok(W,"=");
			token= strtok(NULL,"=");
			netlist[numeroElementos].W = strtod(token, &ptr);

			if (strcmp (netlist[numeroElementos].nome1, "PMOS") == 0)
			{
				netlist[numeroElementos].lambda *= -1;
				netlist[numeroElementos].K *= -1;
			}

			printf("%s %s %s %s %s %s %f %f %f %f %f %f %f %f\n",
																netlist[numeroElementos].nome,na,nb,nc,nd,netlist[numeroElementos].nome1,
																netlist[numeroElementos].L, netlist[numeroElementos].W, netlist[numeroElementos].K,
																netlist[numeroElementos].Vt0, netlist[numeroElementos].lambda, netlist[numeroElementos].gama,
																netlist[numeroElementos].phi, netlist[numeroElementos].Ld);

			netlist[numeroElementos].d=NumerarNo(na);
			netlist[numeroElementos].g=NumerarNo(nb);
			netlist[numeroElementos].s=NumerarNo(nc);
			netlist[numeroElementos].b=NumerarNo(nd);
		}

		else if (tipo=='*')
		{ /* Comentario comeca com "*" */
			printf("Comentario: %s",txt);
			numeroElementos--;
		}
		else if (tipo=='.') //Comando
		{
			if (strcmp(netlist[numeroElementos].nome, ".AC") == 0)
			{
				sscanf(p,"%3s%d%lg%lg", tipoAC, &numeroPontos,&freqInicial, &freqFinal);
				printf ("Tipo AC: %s\nNumero de Pontos: %i\nFreq Inicial: %g\nFreq Final: %g\n", tipoAC, numeroPontos, freqInicial, freqFinal);

			}
			numeroElementos--;
		}
		else
		{
			printf("Elemento desconhecido: %s\n",txt);
			getch();
			exit(1);
		}
	}
	fclose(arquivo);


	/* Acrescenta variaveis de corrente acima dos nos, anotando no netlist */

	// Essa variavel numeroVariaveis é modificada nas chamadas de NumerarNo - Soma 1 a cada novo no
	numeroNos=numeroVariaveis;
	for (i=1; i<=numeroElementos; i++)
	{
		tipo=netlist[i].nome[0];
		if (tipo=='V' || tipo=='E' || tipo=='F' || tipo=='O' || tipo=='L')
		{
			numeroVariaveis++;
			if (numeroVariaveis>MAX_NOS)
			{
				printf("As correntes extra excederam o numero de variaveis permitido (%d)\n",MAX_NOS);
				exit(1);
			}
			strcpy(lista[numeroVariaveis],"j"); /* Tem espaco para mais dois caracteres */
			strcat(lista[numeroVariaveis],netlist[i].nome);
			netlist[i].x=numeroVariaveis;
		}
		else if (tipo=='H')
		{
			numeroVariaveis=numeroVariaveis+2;
			if (numeroVariaveis>MAX_NOS)
			{
				printf("As correntes extra excederam o  numero de variaveis permitido (%d)\n",MAX_NOS);
				exit(1);
			}
			strcpy(lista[numeroVariaveis-1],"jx"); strcat(lista[numeroVariaveis-1],netlist[i].nome);
			netlist[i].x=numeroVariaveis-1;
			strcpy(lista[numeroVariaveis],"jy"); strcat(lista[numeroVariaveis],netlist[i].nome);
			netlist[i].y=numeroVariaveis;
		}
	}
	getch();
	/* Lista tudo */
	printf("Variaveis internas: \n");
	for (i=0; i<=numeroVariaveis; i++)
		printf("%d -> %s\n",i,lista[i]);
	getch();
	printf("Netlist interno final\n");
	for (i=1; i<=numeroElementos; i++) {
		tipo=netlist[i].nome[0];
		if (tipo=='R' || tipo=='L' || tipo=='C')
		{
			printf("%s %d %d %g\n",netlist[i].nome,netlist[i].a,netlist[i].b,netlist[i].valor);
		}
		else if (tipo=='G' || tipo=='E' || tipo=='F' || tipo=='H')
		{
			printf("%s %d %d %d %d %g\n",netlist[i].nome,netlist[i].a,netlist[i].b,netlist[i].c,netlist[i].d,netlist[i].valor);
		}
		else if (tipo=='I' || tipo=='V')
		{
			printf("%s %d %d %g %g %g\n",netlist[i].nome,netlist[i].a,netlist[i].b,netlist[i].modulo,netlist[i].fase,netlist[i].valor);
		}
		else if (tipo=='K')
		{
			printf("%s %s %s %g\n",netlist[i].nome,netlist[i].nome1,netlist[i].nome2,netlist[i].valor);
		}
		else if (tipo=='O')
		{
			printf("%s %d %d %d %d\n",netlist[i].nome,netlist[i].a,netlist[i].b,netlist[i].c,netlist[i].d);
		}
		if (tipo=='V' || tipo=='E' || tipo=='F' || tipo=='O')
			printf("Corrente jx: %d\n",netlist[i].x);
		else if (tipo=='H')
			printf("Correntes jx e jy: %d, %d\n",netlist[i].x,netlist[i].y);
	}
	getch();
	/* Monta o sistema nodal modificado */
	printf("O circuito tem %d nos, %d variaveis e %d elementos\n",numeroNos,numeroVariaveis,numeroElementos);

	getch();

	if (sistemaLinear == 1)
	{
		ZerarSistema();

		/* Monta estampas */
		MontarEstampasPontoOp();

		PrintarSistema(Yn);

		for (i=1; i<=numeroVariaveis; i++)
		{
			for (j=1; j<=numeroVariaveis+1; j++)
				if (YnAnterior[i][j]!=0) printf("%+3.1f ",YnAnterior[i][j]);
				else printf(" ... ");
			printf("\n");
		}

		//return 1;

		/* Resolve o sistema */
		if (ResolverSistema())
		{
			getch();
			exit(1);
		}
	}
	else
	{
		NewtonRaphsonPontoOp();
		return 0;
	}

	#ifdef DEBUG
		 Opcional: Mostra o sistema resolvido
		printf("Sistema resolvido:\n");
		PrintarSistema();
		getch();
	#endif

	//Mostra solucao
	printf("Solucao:\n");
	strcpy(txt,"Tensao");
	for (i=1; i<=numeroVariaveis; i++)
	{
		if (i==numeroNos+1) strcpy(txt,"Corrente");
		printf("%s %s: %g\n",txt,lista[i],Yn[i][numeroVariaveis+1]);
	}
	//getch();

	nomeArquivoTab = strtok(nomeArquivo,".");
	strcat (nomeArquivoTab,".tab");
	arquivoTab = fopen (nomeArquivoTab, "w");
	if (arquivoTab==0)
	{
		printf("Falha ao criar arquivo %s\n",nomeArquivoTab);
		exit (1);
	}
	fprintf (arquivoTab, "f ");
	for (i=1; i<=numeroVariaveis; i++)
	{
		//if (i==numeroNos+1)
		fprintf(arquivoTab, "%sm %sf ",lista[i],lista[i]);
	}
	fprintf (arquivoTab, "\n");

	if (strcmp(tipoAC, "DEC") == 0)
	{
		deltaOmega = pow (10, 1/(double)numeroPontos);

		for (ponto=0; ponto<(int) (numeroPontos*log10(freqFinal/freqInicial))+1; ponto++)
		{


			// Trocar para radianos!!!
			omega = freqInicial * pow (deltaOmega,ponto)*2*M_PI;
			//printf ("Omega = %g\n", omega);

			ZerarSistemaComplexo();
			//PegarPontoOp
			MontarEstampasAnaliseFrequencia(omega);

			if (ResolverSistemaComplexo())
			{
				getch();
				exit(1);
			}

			/*
		#ifdef DEBUG
			//Opcional: Mostra o sistema resolvido
			printf("Sistema resolvido:\n");
			PrintarSistema();
			getch();
		#endif
			 */

			fprintf (arquivoTab, "%g ", omega/(2*M_PI));

			//Mostra solucao
			//printf("Solucao:\n");
			for (i=1; i<=numeroVariaveis; i++)
			{
				fprintf (arquivoTab, "%g %g ", cabs(YnComplex[i][numeroVariaveis+1]), carg(YnComplex[i][numeroVariaveis+1]));
			}
			fprintf(arquivoTab, "\n");
		}
	}
	else if (strcmp(tipoAC,"LIN") == 0)
	{
		deltaOmega = (freqFinal - freqInicial)/(float)numeroPontos;
		deltaOmega = deltaOmega * 2*M_PI;

		for (ponto = 0; ponto <= numeroPontos; ponto++)
		{
			omega = freqInicial*2*M_PI + ponto*deltaOmega;

			ZerarSistemaComplexo();
			//PegarPontoOp
			MontarEstampasAnaliseFrequencia (omega);

			if (ResolverSistemaComplexo())
			{
				getch();
				exit(1);
			}

			fprintf (arquivoTab, "%g ", omega/(2*M_PI));

			for (i=1; i<=numeroVariaveis; i++)
			{
				fprintf (arquivoTab, "%g %g ", cabs(YnComplex[i][numeroVariaveis+1]), carg(YnComplex[i][numeroVariaveis+1]));
			}
			fprintf(arquivoTab, "\n");
		}
	}

	fclose (arquivoTab);

	//getch();

	return 0;
}

