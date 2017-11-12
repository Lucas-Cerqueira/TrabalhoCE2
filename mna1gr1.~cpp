//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "mna1gr1.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TForm1 *Form1;

#include <stdio.h>
#include <conio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include <math.h>
#include <complex.h>
#define MAX_LINHA 200
#define MAX_NOME 11
#define MAX_ELEM 50
#define MAX_NOS 50
#define TOLG 1e-12
//#define DEBUG

#define RESISTENCIA 1e9

#define MAX_ITERACOES 500
#define MAX_INIT 100
#define ERRO_MAX 1e-6
#define X_ERRO 1

#define I complex <double> (0.0,1.0)

using namespace std;

#define printf xprintf
void xprintf(char* format,...) /* Escreve no memo1 */
{
  va_list paramlist;
  char txt[2000];
  va_start(paramlist,format);
  vsprintf(txt,format,paramlist);
  Form1->Memo1->Lines->Add(txt);
}

typedef enum TipoModoOperacao
{
	corte, triodo, saturacao
} TipoModoOperacao;

typedef struct elemento // Elemento do netlist
{
	char nome[MAX_NOME], nome1[MAX_NOME], nome2[MAX_NOME];
	double valor, modulo, fase;
	double L, W, K, Vt0, lambda, gama, phi, Ld, gm, gds, gmb, i0, cgs, cgd, cgb;
	int a,b,c,d,x,y,g,s;
	TipoModoOperacao modoOperacao;
} elemento;


elemento netlist[MAX_ELEM]; /* Netlist */

bool
mostrarEstampasDC,
mostrarEstampasAC,
mostrarParametrosMOS,
iniciarConduzindo;

int
numeroElementos, /* Numero de elementos */
numeroVariaveis, /* Numero de variaveis */
numeroNos, /* Numero de nos */
i,j,k,x,y,
ponto,
numeroPontos,
sistemaLinear,
analiseFreq,
iteracao, init;


char
/* Foram colocados limites nos formatos de leitura para alguma protecao
   contra excesso de caracteres nestas variaveis */
nomeArquivo[MAX_LINHA+1],
*nomeArquivoTab,
mosResult[MAX_LINHA],
mosGm[MAX_LINHA],
tipo,
na[MAX_NOME],nb[MAX_NOME],nc[MAX_NOME],nd[MAX_NOME],
lista[MAX_NOS+1][MAX_NOME+2],
txt[MAX_LINHA+1],
*p,
tipoAC[4];

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

complex <double>
gComplex,
YnComplex[MAX_NOS+1][MAX_NOS+2];

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
	g = 1/RESISTENCIA;

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

int MontarCapacitanciasMOSAnaliseFreq(elemento* mos, double w)
{
        complex <double> g;
        
	//Capacitancia G-S
	g = I*w*mos->cgs;
	YnComplex[mos->s][mos->s] += g;
	YnComplex[mos->g][mos->g] += g;
	YnComplex[mos->s][mos->g] -= g;
	YnComplex[mos->g][mos->s] -= g;

	//Capacitancia G-D
	g = I*w*mos->cgd;
	YnComplex[mos->d][mos->d] += g;
	YnComplex[mos->g][mos->g] += g;
	YnComplex[mos->d][mos->g] -= g;
	YnComplex[mos->g][mos->d] -= g;

	//Capacitancia G-B
	g = I*w*mos->cgb;
	YnComplex[mos->b][mos->b] += g;
	YnComplex[mos->g][mos->g] += g;
	YnComplex[mos->b][mos->g] -= g;
	YnComplex[mos->g][mos->b] -= g;

	return 0;
}

int MontarEstampaMOSPontoOp (elemento* mos)
{
	double vgs,vds,vbs,Vt, Cox, mi;
	int auxiliar;


	if (mos->nome1[0]=='N')
		mi = 0.05;
	else
		mi = 0.02;

	Cox = 2*mos->K/mi;

	// Garante que D é o nó de maior tensão no NMOS
	if (mos->nome1[0]=='N')
	{
		if (YnAnterior[mos->d][numeroVariaveis+1] < YnAnterior[mos->s][numeroVariaveis+1])
		{
			auxiliar = mos->d;
                        mos->d = mos->s;
			mos->s = auxiliar;
		}
	}
	// Garante que S é o nó de maior tensão no PMOS
	else
	{
		if (YnAnterior[mos->d][numeroVariaveis+1] > YnAnterior[mos->s][numeroVariaveis+1])
		{
			auxiliar = mos->d;
                        mos->d = mos->s;
			mos->s = auxiliar;
		}
	}


	vgs = YnAnterior[mos->g][numeroVariaveis+1] - YnAnterior[mos->s][numeroVariaveis+1];
	vds = YnAnterior[mos->d][numeroVariaveis+1] - YnAnterior[mos->s][numeroVariaveis+1];
	vbs = YnAnterior[mos->b][numeroVariaveis+1] - YnAnterior[mos->s][numeroVariaveis+1];

        // Limitando Vbs
        if ((mos->nome1[0]=='N') && (vbs > 0))
        	vbs = 0;
        else if ((mos->nome1[0]=='P') && (-vbs > 0))
                vbs = 0;

        // Cálculo de Vt levando em conta efeito de Vbs
	Vt = mos->Vt0 + mos->gama*(sqrt(mos->phi - vbs) - sqrt(mos->phi));

        // Inicia o transistor fora do corte se a opção for marcada
	if (iteracao == 1 && init == 1 && iniciarConduzindo)
	{
		vgs = 2*Vt;
	}

        if (mos->nome1[0]=='P')
	{
		vgs *= -1;
		vds *= -1;
                vbs *= -1;
        }

        if (vgs<Vt) // Corte
        {
                mos->modoOperacao = corte;

                mos->gds = mos->gm = mos->gmb =mos->i0 = 0;

                mos->cgd = mos->cgs = Cox*mos->W*mos->Ld;
                mos->cgb = Cox*mos->W*mos->L;
        }
        else if (vds <= (vgs-Vt)) // Triodo
        {
                mos->modoOperacao = triodo;

                mos->gm = mos->K * (mos->W/mos->L) * (2*vds) * (1 + mos->lambda*vds);
                mos->gds = mos->K * (mos->W/mos->L) * (2*(vgs-Vt) - 2*vds + 4*mos->lambda*(vgs-Vt)*vds - 3*mos->lambda*vds*vds);
                mos->gmb = (mos->gm*mos->gama)/(2*sqrt(mos->phi-vbs));

                mos->i0 = mos->K * (mos->W/mos->L) * (2*(vgs-Vt)*vds - vds*vds) * (1+mos->lambda*vds) - mos->gm*vgs - mos->gds*vds - mos->gmb*vbs;

                mos->cgs = Cox * mos->W * mos->Ld + Cox * mos->W * mos->L/2.0;
                mos->cgd = mos->cgs;
                mos->cgb = 0;
        }
        else // Saturacao
        {
                mos->modoOperacao = saturacao;

                mos->gm = mos->K * (mos->W/mos->L) * 2 * (vgs-Vt) * (1+mos->lambda*vds);
                mos->gds = mos->K * (mos->W/mos->L) * mos->lambda * (vgs-Vt)*(vgs-Vt);
                mos->gmb = (mos->gm*mos->gama)/(2*sqrt(mos->phi-vbs));

                mos->i0 = mos->K * (mos->W/mos->L) * (vgs-Vt)*(vgs-Vt) * (1+mos->lambda*vds) - mos->gm*vgs - mos->gds*vds - mos->gmb*vbs;

                mos->cgs = Cox * mos->W * mos->Ld + 2*Cox * mos->W * mos->L/3.0;
                mos->cgd = Cox * mos->W * mos->Ld;
                mos->cgb = 0;
        }

	if (mos->nome1[0]=='P')
	{
		mos->i0 *= -1;
	}

	Yn[mos->d][mos->d] += mos->gds;
	Yn[mos->d][mos->s] -= mos->gds + mos->gm + mos->gmb;
	Yn[mos->d][mos->g] += mos->gm;
	Yn[mos->d][mos->b] += mos->gmb;
	Yn[mos->s][mos->d] -= mos->gds;
	Yn[mos->s][mos->s] += mos->gds + mos->gm + mos->gmb;
	Yn[mos->s][mos->g] -= mos->gm;
	Yn[mos->s][mos->b] -= mos->gmb;

	Yn[mos->d][numeroVariaveis+1] -= mos->i0;
	Yn[mos->s][numeroVariaveis+1] += mos->i0;



        if (mostrarParametrosMOS)
        {
                if (mos->nome1[0]=='P')
        	{
	        	vgs *= -1;
        		vds *= -1;
                        vbs *= -1;
                }
                printf ("\nn=%i | %s Vd:%g Vg:%g Vs:%g Vb:%g Vgs:%g Vds:%g Vbs:%g Vt:%g Gmb:%g Gm:%g Gds:%g i0:%g\n", iteracao, mos->nome,
                        YnAnterior[mos->d][numeroVariaveis+1], YnAnterior[mos->g][numeroVariaveis+1], YnAnterior[mos->s][numeroVariaveis+1],
                        YnAnterior[mos->b][numeroVariaveis+1], vgs, vds, vbs, Vt, mos->gmb, mos->gm, mos->gds, mos->i0);
        }
	return 0;
}

void PrintarSistema(double sistema[MAX_NOS+1][MAX_NOS+2])
{
        int i,j;
        char linha[2000], resultado[20];
	for (i=1; i<=numeroVariaveis; i++)
	{
                strcpy (linha, "");
		for (j=1; j<=numeroVariaveis+1; j++)
                {
			if (sistema[i][j]!=0)
                        {
                                sprintf (resultado,"%+3.2e ",sistema[i][j]);
                                strcat(linha, resultado);
                        }
			else strcat(linha,"......... ");
                }
		printf(linha);
	}
}

void PrintarSistemaComplexo()
{
        int k,j;
        char linha[2000], resultado[20];
	for (k=1; k<=numeroVariaveis; k++)
	{
                strcpy (linha, "");
		for (j=1; j<=numeroVariaveis+1; j++)
                {
			if (abs(YnComplex[k][j]) != 0)
                        {
                                sprintf(resultado,"%+3.2f%+3.2fj ", real(YnComplex[k][j]), imag(YnComplex[k][j]));
				strcat(linha,resultado);
                        }
			else strcat(linha," ........... ");
                }
		printf(linha);
	}
        printf ("");
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
		if (fabs(t)<TOLG)
		{
			printf("Sistema singular\n");
			PrintarSistema(Yn);
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
	complex <double> t, p;

	for (i=1; i<=numeroVariaveis; i++) {
		t=0.0;
		a=i;
		for (l=i; l<=numeroVariaveis; l++) {
			if (abs(YnComplex[l][i])>abs(t)) {
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
		if (abs(t)<TOLG)
		{
			printf("Sistema singular\n");
			PrintarSistemaComplexo();
			return 1;
		}
		for (j=numeroVariaveis+1; j>0; j--) {  /* Basta j>i em vez de j>0 */
			YnComplex[i][j]/= t;
			p=YnComplex[i][j];
			if (p!=(0.0,0.0))  /* Evita operacoes com zero */
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
	while (!achou && i<=numeroVariaveis)
		if (!(achou=!strcmp(nome,lista[i])))
			i++;

	if (!achou)
	{
		if (numeroVariaveis==MAX_NOS)
		{
			printf("O programa so aceita ate %d nos\n",numeroVariaveis);
			return -1;
		}
		numeroVariaveis++;
		strcpy(lista[numeroVariaveis],nome);
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
			g = 1/1e-9;
			Yn[netlist[i].a][netlist[i].x] += 1;
			Yn[netlist[i].b][netlist[i].x] -= 1;
			Yn[netlist[i].x][netlist[i].a] -= 1;
			Yn[netlist[i].x][netlist[i].b] += 1;
			Yn[netlist[i].x][netlist[i].x] += 1/g;
		}
		else if (tipo=='C')
		{
			g = 1/(double)RESISTENCIA;
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
		else if (tipo=='I')  // I e V sÃ³ o valor contÃ­nuo
		{
			g = netlist[i].valor;
			Yn[netlist[i].a][numeroVariaveis+1] -= g;
			Yn[netlist[i].b][numeroVariaveis+1] += g;
		}
		else if (tipo=='V')
		{
			Yn[netlist[i].a][netlist[i].x] += 1;
			Yn[netlist[i].b][netlist[i].x] -= 1;
			Yn[netlist[i].x][netlist[i].a] += 1;
			Yn[netlist[i].x][netlist[i].b] -= 1;
			Yn[netlist[i].x][numeroVariaveis+1] += netlist[i].valor;
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
			MontarEstampaMOSPontoOp(&netlist[i]);
			MontarCapacitanciasMOSPontoOp(netlist[i]);
		}
                if (mostrarEstampasDC)
                {
        		printf("Sistema após a estampa de %s (n=%i):",netlist[i].nome,iteracao);
        		PrintarSistema(Yn);
                }
	}
}

void MontarEstampasAnaliseFrequencia(double w)
{
	for (i=1; i<=numeroElementos; i++)
	{
		tipo=netlist[i].nome[0];

		if (tipo=='M')
		{
			YnComplex[netlist[i].d][netlist[i].d] += netlist[i].gds;
			YnComplex[netlist[i].d][netlist[i].s] -= netlist[i].gds + netlist[i].gm + netlist[i].gmb;
			YnComplex[netlist[i].d][netlist[i].g] += netlist[i].gm;
			YnComplex[netlist[i].d][netlist[i].b] += netlist[i].gmb;
			YnComplex[netlist[i].s][netlist[i].d] -= netlist[i].gds;
			YnComplex[netlist[i].s][netlist[i].s] += netlist[i].gds + netlist[i].gm + netlist[i].gmb;
			YnComplex[netlist[i].s][netlist[i].g] -= netlist[i].gm;
			YnComplex[netlist[i].s][netlist[i].b] -= netlist[i].gmb;

			MontarCapacitanciasMOSAnaliseFreq (&netlist[i], w);
		}
		else if (tipo=='R')
		{
			gComplex = 1/netlist[i].valor;
			YnComplex[netlist[i].a][netlist[i].a] += gComplex;
			YnComplex[netlist[i].b][netlist[i].b] += gComplex;
			YnComplex[netlist[i].a][netlist[i].b] -= gComplex;
			YnComplex[netlist[i].b][netlist[i].a] -= gComplex;
		}
		else if (tipo=='L')
		{
			gComplex = 1./(I*w*netlist[i].valor); // 1/jwL
			//printf ("Modulo gL = %g\n", cabs (gComplex));
			YnComplex[netlist[i].a][netlist[i].x] += 1;
			YnComplex[netlist[i].b][netlist[i].x] -= 1;
			YnComplex[netlist[i].x][netlist[i].a] -= 1;
			YnComplex[netlist[i].x][netlist[i].b] += 1;
			YnComplex[netlist[i].x][netlist[i].x] += 1./gComplex;
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
			//printf ("Capacitor: %g + %gi\n", real(gComplex), imag(gComplex));
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
                if (mostrarEstampasAC)
                {
        		printf("Sistema apos a estampa de %s (omega=%f)\n",netlist[i].nome,w);
	        	PrintarSistemaComplexo();
                }
	}
}

int TestarConvergencia()
{
	int i, convergiu;

        convergiu = 1;

	for (i=1; i<=numeroVariaveis;i++)
	{
		if (fabs(Yn[i][numeroVariaveis+1]) > X_ERRO) // Usa erro relativo
		{
			erros[i] = X_ERRO*fabs((Yn[i][numeroVariaveis+1]-YnAnterior[i][numeroVariaveis+1])/Yn[i][numeroVariaveis+1]);
		}
		else // Usa erro absoluto
		{
			erros[i] = fabs(Yn[i][numeroVariaveis+1]-YnAnterior[i][numeroVariaveis+1]);
		}

                /*if (i <= numeroNos)
                        printf ("Tensao %i: %4.3e   |  Erro: %4.3e", i, Yn[i][numeroVariaveis+1], erros[i]);
                else
                        printf ("Corrente %i: %4.3e |  Erro: %4.3e", i, Yn[i][numeroVariaveis+1], erros[i]);*/
                

		if (erros[i] > ERRO_MAX)
			convergiu = 0;
	}
	return convergiu;
}

void InicializacaoRandomica()
{
	double valor;

	srand ((unsigned)time(NULL));

	for (i=1; i<=numeroVariaveis;i++)
	{
		if (erros[i] > ERRO_MAX)
		{
			if (i <= numeroNos)  // Se for tensao, randomiza de -10 a 10
			{
				valor = rand() % 20001;
				valor -= 10000;
				valor /= 1000.0;
			}
			else                // Se for corrente, randomiza de -1 a 1
			{
				valor = rand() % 2001;
				valor -= 1000;
				valor /= 1000.0;
			}
			YnAnterior[i][numeroVariaveis+1] = valor;
		}
		else    //Mantem a solucao encontrada
			YnAnterior[i][numeroVariaveis+1] = Yn[i][numeroVariaveis+1];
	}
}

int NewtonRaphsonPontoOp ()
{
	int i, convergiu;

	// Zerando Yn
	ZerarSistema();
	SalvarUltimaIteracao(); // Zerar o YnAnterior tambem

	//Inicializacao das variaveis do sistema
	for (i=1; i<=numeroVariaveis; i++)
	{
		YnAnterior[i][numeroVariaveis+1] = 0.1;
	}

	for (init=1; init<=MAX_INIT; init++)
	{
                printf ("Inicializacao: %i", init);
		if (init != 1)
		{
                        InicializacaoRandomica();
		}

		for (iteracao = 0; iteracao<=MAX_ITERACOES; iteracao++)
		{
			ZerarSistema();
			MontarEstampasPontoOp();
			if (ResolverSistema())
                                return 2;

                        //PrintarSistema(Yn);

			convergiu = TestarConvergencia();
                        SalvarUltimaIteracao();
			if (convergiu == 1)
			{
				printf ("Sistema convergiu com %i iteracoes e %i inicializacoes randomicas\n",iteracao, init-1);
				return 0;
			}
		}
	}
        printf ("Não convergiu");
        return 1;
}

int AnalisarFrequencia (double deltaOmega, bool linear)
{
        omega = freqInicial*2*M_PI;
        while (omega <= (freqFinal*2*M_PI + 0.000001))
        {
                ZerarSistemaComplexo();
		MontarEstampasAnaliseFrequencia(omega);
                if (ResolverSistemaComplexo())
			return 1;

                fprintf (arquivoTab, "%g ", omega/(2*M_PI));
		for (i=1; i<=numeroVariaveis; i++)
		{
                        if (abs(YnComplex[i][numeroVariaveis+1])!= 0)
                                fprintf (arquivoTab, "%g %g", abs(YnComplex[i][numeroVariaveis+1]), arg(YnComplex[i][numeroVariaveis+1])*180/M_PI);

                        else
                                fprintf (arquivoTab, "0 0");

                        if (i!=numeroVariaveis) // Se não for a ultima, coloca espaço entre as variaveis
                                fprintf (arquivoTab, " ");
		}
		fprintf(arquivoTab, "\n");

                if (linear)
                        omega += deltaOmega;
                else
        		omega *= deltaOmega;

        }
        
       	fclose (arquivoTab);
        return 0;
}

//---------------------------------------------------------------------------
__fastcall TForm1::TForm1(TComponent* Owner)
        : TForm(Owner)
{
}
//---------------------------------------------------------------------------
void __fastcall TForm1::Abrir1Click(TObject *Sender)
{
        denovo:
        
        if (MostrarEstampasDC1->Checked)
                mostrarEstampasDC = true;
        else
                mostrarEstampasDC = false;

        if (MostrarEstampasAC1->Checked)
                mostrarEstampasAC = true;
        else
                mostrarEstampasAC = false;

        if (IniciarMOSconduzindo1->Checked)
                iniciarConduzindo = true;
        else
                iniciarConduzindo = false;

        if (MostrarParametrosMOS1->Checked)
                mostrarParametrosMOS = true;
        else
                mostrarParametrosMOS = false;

        analiseFreq = 0;
	sistemaLinear = 1;
	numeroElementos=0;
	numeroVariaveis=0;
	strcpy(lista[0],"0");
  /* Leitura do netlist */
        if (!OpenDialog1->Execute()) return;
        strcpy(nomeArquivo,OpenDialog1->FileName.c_str());
        arquivo=fopen(nomeArquivo,"r");
        if (arquivo==0)
        {
            printf("Arquivo %s inexistente",nomeArquivo);
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
			goto fim;
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
                        netlist[numeroElementos].modulo = netlist[numeroElementos].fase = netlist[numeroElementos].valor = 0;
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

                        // Cálculo do M
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
			//*Transistor MOS: M<nome> <nÃ³d> <nÃ³g> <nÃ³s> <nÃ³b> <NMOS ou PMOS> L=<comprimento> W=<largura> <K> <Vt0> <lambda> <gama> <phi> <Ld>
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
			token = strtok(NULL,"=");
			netlist[numeroElementos].W = strtod(token, &ptr);


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
				analiseFreq = 1;

				sscanf(p,"%3s%d%lg%lg", tipoAC, &numeroPontos,&freqInicial, &freqFinal);
				printf ("Tipo AC: %s\nNumero de Pontos: %i\nFreq Inicial: %g\nFreq Final: %g\n", tipoAC, numeroPontos, freqInicial, freqFinal);

			}
			numeroElementos--;
		}
		else
		{
			printf("Elemento desconhecido: %s\n",txt);
			getch();
			goto fim;
		}
	}
	fclose(arquivo);


	/* Acrescenta variaveis de corrente acima dos nos, anotando no netlist */

	// Essa variavel numeroVariaveis Ã© modificada nas chamadas de NumerarNo - Soma 1 a cada novo no
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
				goto fim;
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
				goto fim;
			}
			strcpy(lista[numeroVariaveis-1],"jx"); strcat(lista[numeroVariaveis-1],netlist[i].nome);
			netlist[i].x=numeroVariaveis-1;
			strcpy(lista[numeroVariaveis],"jy"); strcat(lista[numeroVariaveis],netlist[i].nome);
			netlist[i].y=numeroVariaveis;
		}
	}

	/* Monta o sistema nodal modificado */
	printf("O circuito tem %d nos, %d variaveis e %d elementos\n",numeroNos,numeroVariaveis,numeroElementos);

	if (sistemaLinear == 1)
	{
		ZerarSistema();

		MontarEstampasPontoOp();

		if (ResolverSistema())
			goto fim;

	}
	else if (NewtonRaphsonPontoOp() != 0)
                goto fim; //Problema na convergencia


	//Mostra solucao
	printf("Solucao:\n");
	strcpy(txt,"Tensao");
	for (i=1; i<=numeroVariaveis; i++)
	{
		if (i==numeroNos+1) strcpy(txt,"Corrente");
		printf("%s %s: %g\n",txt,lista[i],Yn[i][numeroVariaveis+1]);
	}
	for (i=1; i<=numeroElementos; i++)
	{
		if (netlist[i].nome[0] == 'M')
		{
        	        sprintf (mosResult,"%s: %s ", netlist[i].nome, netlist[i].nome1);

			if (netlist[i].modoOperacao == corte)
				strcat (mosResult,"cortado ");
			else if (netlist[i].modoOperacao == triodo)
				strcat (mosResult,"triodo ");
			else
				strcat (mosResult,"saturado ");

			double vgs,vds,vbs;
                        // Pega os últimos resultados
                        vgs = YnAnterior[netlist[i].g][numeroVariaveis+1] - YnAnterior[netlist[i].s][numeroVariaveis+1];
			vds = YnAnterior[netlist[i].d][numeroVariaveis+1] - YnAnterior[netlist[i].s][numeroVariaveis+1];
			vbs = YnAnterior[netlist[i].b][numeroVariaveis+1] - YnAnterior[netlist[i].s][numeroVariaveis+1];

			sprintf (mosGm,"Gm=%e Gds=%e Gmb=%e Cgs=%e Cgd=%e Cgb=%f Vgs=%g Vds=%g Vbs=%g Id=%g",
                                        netlist[i].gm, netlist[i].gds, netlist[i].gmb, netlist[i].cgs, netlist[i].cgd, netlist[i].cgb,
                                        vgs, vds, vbs, netlist[i].i0 + netlist[i].gm*vgs + netlist[i].gds*vds + netlist[i].gmb*vbs);
                        strcat (mosResult,mosGm);
                        printf (mosResult);
		}
	}

	if (analiseFreq == 1)
	{
		printf ("");
		printf ("Iniciando Analise em frequencia");

		nomeArquivoTab = strtok(nomeArquivo,".");
		strcat (nomeArquivoTab,".tab");
		arquivoTab = fopen (nomeArquivoTab, "w");
		if (arquivoTab==0)
		{
			printf("Falha ao criar arquivo %s",nomeArquivoTab);
			goto fim;
		}
		fprintf (arquivoTab, "f ");
		for (i=1; i<=numeroVariaveis; i++)
		{
			fprintf(arquivoTab, "%sm %sf ",lista[i],lista[i]);
		}
		fprintf (arquivoTab, "\n");

		if (strcmp(tipoAC, "DEC") == 0)
		{
			deltaOmega = std::pow (10, 1/(double)(numeroPontos-1));
                        AnalisarFrequencia (deltaOmega, false);

		}
		else if (strcmp(tipoAC, "OCT") == 0)
		{
			deltaOmega = std::pow (2, 1/(double)(numeroPontos-1));
                        AnalisarFrequencia (deltaOmega, false);
		}
		else if (strcmp(tipoAC,"LIN") == 0)
		{
			deltaOmega = (freqFinal - freqInicial)/(float)numeroPontos;
			deltaOmega = deltaOmega * 2*M_PI;

                        AnalisarFrequencia (deltaOmega, true);
		}
	}

        printf ("Resultado escrito no arquivo %s", nomeArquivoTab);

        fim:
}

