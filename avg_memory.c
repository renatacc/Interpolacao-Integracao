/* 
 * Trabalho de Cálculo Numérico
 * Author: Ana Paula Fernandes de Souza - 0022647
 *         Renata Caroline Cunha - 0022429
 * Data: 05/07/2016
 * 
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
//Estrutura criada para armazenar os valores dos pontos
struct ponto{
                double x,y;    
};
struct dados{
                int num_pares; //Número de pontos
                double max_x,min_x; //Valores máximo e mínimo de x
                double max_y,min_y;//Valores máximo e mínimo de y
                struct ponto *pares; //Ponteiro para os pontos
};  

typedef struct dados *Pontos;  //Ponteiro para estrutura

//Requisito 01 - TAD para Sequência Pares Ordenados
Pontos CriaEstrutura(char *nome) 
{
    FILE *fp;
    int lixo,tam;
    fp=fopen(nome,"r");

    if(fp==NULL)
    {
        return(0);
    }
    while(fscanf(fp,"%d %d",&tam,&lixo)!=EOF)
    {
        //não vai fazer nada.
    }    
    fclose(fp);
    
    tam++;
    Pontos aux=(struct dados*)malloc(sizeof(struct dados));  //Alocação da memória
    aux->num_pares=tam; 
    aux->pares=(struct ponto*)malloc((tam+1)*sizeof(struct ponto)); //Estrutura indexada de 1, por isso tam+1.
    aux->max_x=0;
    aux->min_x=0;
    aux->max_y=0;
    aux->min_y=0;

    return(aux);   
}

//Requisito 02 - Entrada de dados
void CarregaDados(char* nome,Pontos p)
{   
  FILE *fp;
  double x,y;
  int i;
  fp=fopen(nome,"r");

  if(fp==NULL)
  {
      return;
  }
  for(i=1;i<=p->num_pares;i++)
  {
        fscanf(fp,"%lf %lf",&x,&y);  //Lendo os dados do arquivo e armazenando no vetor pontos
        p->pares[i].x=x;  
        p->pares[i].y=y;      
        
        //Verificando os máximos e mínimos
        if(i==1)  //Caso base
        {
            p->min_x=x;
            p->max_x=x;
            
            p->min_y=y;
            p->max_y=y;
        }
        else //Outros casos
        {
            if(x<p->min_x)
            {
                 p->min_x=x;
            }
            else
            {
                if(x>p->max_x)
                {
                    p->max_x=x;
                }      
            }
    
            if(y<p->min_y)
            {
                p->min_y=y;
            }
            else
            {
                if(y>p->max_y)
                {
                    p->max_y=y;
                }            
            }        
        }
  }  
  fclose(fp); 
}
    
//Requisito 03 - Cálculo das 2ª Derivada da Spline Cúbica
double* CalculaDerivadaSpline(Pontos p) 
{
    if(p->num_pares<3)
    {
        return(0); //Como não possui 3 pontos retornaremos o valor 0.
    }
    
    //Sistema Tridiagonal Simétrico
    int m=p->num_pares-2;
    
    double Ha=(p->pares[2].x)-(p->pares[1].x);   
    double DeltaA=(p->pares[2].y-p->pares[1].y)/Ha; 
    
    int i;
    
    double *e=(double*)malloc((m+1)*sizeof(double));
    double *d=(double*)malloc((m+1)*sizeof(double));
    double *s2=(double*)malloc((p->num_pares+1)*sizeof(double));
    
    for(i=1; i<=m; i++)
    {
         double Hb=(p->pares[i+2].x)-(p->pares[i+1].x);
         double DeltaB=(p->pares[i+2].y-p->pares[i+1].y)/Hb;
         
         e[i]=Hb;
         d[i]=2*(Ha+Hb);        
         s2[i+1]=6*(DeltaB-DeltaA);
         
         Ha=Hb;
         
         DeltaA=DeltaB;        
    }   
    
    double t;
    //Eliminação de Gauss
    for(i=2; i<=m; i++)
    {
        t=e[i-1]/d[i-1];
        d[i]=d[i]-t*e[i-1];
        s2[i+1]=s2[i+1]-t*s2[i];
    }
    
    //Substituição Retroativa
    s2[m+1]=s2[m+1]/d[m];
    
    for(i=m; i>=2; i--)
    {
        s2[i]=(s2[i]-e[i-1]*s2[i+1])/d[i-1];
    }
    
    s2[1]=0;
    s2[m+2]=0;
    
    free(e); //Liberando a memória
    free(d);
    
    return(s2); 
}

//Requisito 04 - Função que Avalia P(x) usando Splines
double AvaliaSpline(Pontos p, double* s2, double valor)
{
    if(valor<p->pares[1].x || valor>p->pares[p->num_pares].x)
    {
        return(0);
    }   
    //Busca Binária
    int inf=1;
    int sup=p->num_pares;
    
    while((sup-inf)>1)
    {
        int indice=(inf+sup)/2;
        
        if(p->pares[indice].x>valor)
        {
            sup=indice;
        }
        else
        {
            inf=indice;
        }
    }
    
    //Avaliação de Ṕ(x) por Horner
    double h=(p->pares[sup].x)-(p->pares[inf].x);
    double a=(s2[sup]-s2[inf])/(6*h);
    double b=s2[inf]*0.5;
    double c=(((p->pares[sup].y)-(p->pares[inf].y))/h)-((s2[sup]+2*s2[inf])*h/6);
    double d=p->pares[inf].y;
  
    h=valor-(p->pares[inf].x);
    
    double resultado=((a*h+b)*h+c)*h+d;
    
    return(resultado);
}
//Requisito 05 - Gerador de Números Uniforme U(min,max)
double U(int min, int max)
{   
    int x=rand()%(max+1); //Para que seja possível sortear o valor máximo foi colocado o +1.   
    
    double y=(double)x/max;
    
    y=y*(max-min)+min;
    
    return(y);
}

//Requisito 06 - Integral Numérica por Monte Carlo
double IntegralMonteCarlo(int n,Pontos p,double *s2) //s2 é o vetor gerado no Requisito 03 e n é o numero de pontos a ser gerado
{   
    int i;
    int num_abaixo=0;
    p->min_y=0; 
     srand(time(NULL));

    for(i=1;i<=n;i++) 
    {
        double x=U(p->min_x,p->max_x);
        double y=U(p->min_y,(p->max_y)*1.2); //Aumentando a altura em 20% 
        
        double y1=AvaliaSpline(p,s2,x);
        
        if(y<=y1) 
        {
            num_abaixo++;
        }       
    }
     
    int AreaTotal=((p->max_x)-(p->min_x))*((p->max_y)*1.2)-(p->min_y);
   
    double Area=AreaTotal*((double)num_abaixo/n);
 
    return(Area);
}

//Requisito 07 - Uso Médio de Memória usando TVMI
double TVMI(Pontos p, double *s2, int n)
{
      double MonteCarlo=IntegralMonteCarlo(n,p,s2);
      double x=(1/(double)(p->max_x-p->min_x)); 
      double med=x*MonteCarlo;
      
      return(med);
}

//Requesito 08 - Saída de Dados(Terminal)
void SaidaDados(Pontos p, double *s2, int n)
{    
    double resultado=TVMI(p,s2,n);
    
    printf("Number of Samples: %lf \n",p->max_x+1);
    printf("Average Memory Usage : %lf \n",resultado); 
    printf("Run Rscript saida.r to generate Average Memory Usage Chart");
    
    return;
}

//Requesito 09 - Saída de Dados(R Script)
void SalvaDados(Pontos p, double *s2, int n,char *nome)
{   
   SaidaDados(p,s2,n);
   double resultado=TVMI(p,s2,n);
   
      FILE *arq;
      arq=fopen(nome,"w");
 
    if(arq==NULL)
    {
         //Apenas para evitar erros       
    }
    else
    {
        fprintf(arq,"#");
        fprintf(arq,"Generated automatically by \"avg-memory\" application\n");
        fprintf(arq,"\n\n");

        int i;

        fprintf(arq,"# Original points (x coordinates)\n");
        fprintf(arq,"xorig <- c(\n");
        for(i=1;i<=p->num_pares;i++)
        {
            if(i==p->num_pares)
            {
                fprintf(arq,"%.2lf\n",p->pares[i].x);
            }
            else
            {
                fprintf(arq,"%.2lf,\n",p->pares[i].x);
            }    
        }

        fprintf(arq,");\n\n");
        fprintf(arq,"# Original points (y coordinates)\n");
        fprintf(arq,"yorig <- c(\n");

        for(i=1;i<=p->num_pares;i++)
        {
            if(i==p->num_pares)
            {
                fprintf(arq,"%.2lf\n",p->pares[i].y);
            }
            else
            {
                fprintf(arq,"%.2lf,\n",p->pares[i].y);
            }    
        }    

       fprintf(arq,");\n\n");         
    }  
    fprintf(arq,"# Spline points (x coordinates, sampling interval = 0.01)\n");
    fprintf(arq,"xspl <- c(\n");
      
      double x=p->min_x; 
      while(x<(p->max_x)-0.01)
      {
         fprintf(arq,"%.2lf,\n",x);  
         x=x+0.01;                       
      }
       fprintf(arq,"%.2lf \n",x);
       fprintf(arq,");\n\n"); 
       
    //pontos no y                
     fprintf(arq,"# Spline points (y coordinates, sampling interval = 0.01)\n");
     fprintf(arq,"yspl <- c(\n");
     
      x=p->min_x; 
      while(x<(p->max_x)-0.01)
      {
         double aval=AvaliaSpline(p,s2,x); 
         fprintf(arq,"%lf,\n",aval);  
         x=x+0.01;                       
      }
      
       double aval=AvaliaSpline(p,s2,x);      
       fprintf(arq,"%lf \n",aval);
       fprintf(arq,");\n\n");         
                
  fprintf(arq,"# Average Memory Usage \n");   
  fprintf(arq,"AvgMemory <- %lf ; \n",resultado);    
  fprintf(arq,"# Plot the values in .png file \n");
  fprintf(arq,"png(file=\"saida.png\", width=1200);\n");
  fprintf(arq,"title <- paste(\"AVG Memory Usage: %lf Kb \");\n",resultado);
  fprintf(arq,"plot(xspl, yspl, type=\"l\", col=\"blue\", main=title, xlab=\"Samples\", ylab=\"Mem. Usage\", lwd=3);\n");
  fprintf(arq,"points(xorig, yorig, pch=19, col=\"red\");\n");
  fprintf(arq,"lines(c(min(xorig), max(xorig)), c(AvgMemory, AvgMemory), col=\"black\", lty=2, lwd=3);\n");
  fprintf(arq,"dev.off();\n");
   
 fclose(arq);   
}

int main(int argc, char **argv) 
{
    if(argc!=3)
    {
        return(0);
    }
     
  Pontos p=CriaEstrutura(argv[1]);
          CarregaDados(argv[1],p);

    double *s2=CalculaDerivadaSpline(p);
     
     SalvaDados(p,s2,10000000,argv[2]); //10000000 representa a quantidade de pontos para o monte carlo
  
    free(s2); 
    return (0);

}

