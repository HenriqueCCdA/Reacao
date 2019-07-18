#ifndef _DEFINE_
  #define _DEFINE_
/*...*/
  #include<stdlib.h>
  #include<stdio.h>
/*...................................................................*/

/*...*/
  #define MAX_SPECIES   10
  #define MAX_REAC      10
  #define MAX_NAME_SP   10
  #define R              8.3144621     /*J/mol K*/ 
  #define RC             1.987          /*cal/mol K*/
  #define PREREF         1.01325e+05   /*Pa*/

/*... Kg/kmol*/
  #define MW_O 15.9994e0
  #define MW_H  1.00794e0
  #define MW_N 14.00674e0 
  #define MW_C 12.01100e0
//#define MW_O 16.00000e0
//#define MW_H  1.00000e0
//#define MW_N 14.00000e0 
//#define MW_C 12.00000e0
/*...................................................................*/

/*numero maximo de propriedade*/
  #define MAXPROP      5  /*numero maximo de propriedades*/
  #define MAXMAT      200 /*numero maximo de materias*/
  #define MAX_TRANS_EQ 3 /*numero maximo de equacoes de transporte*/ 
  #define MAX_DIF_EQ   3 /*numero maximo de equacoes de difusa*/ 
/*...................................................................*/

/*...cellLoop*/
  #define  MAX_NUM_RES   28
  #define  MAX_NUM_PONT  168
  #define  MAX_NUM_FACE  6
  #define  MAX_SN        24 
/*...................................................................*/

/*...*/
  #define MAX_NDF 5
/*...................................................................*/

/*... solver*/
  #define PCG        1
  #define PBICGSTAB  2
/*...................................................................*/

/*... CSR*/
  #define CSR  1
  #define CSRD 2
  #define CSRC 3
/*...................................................................*/

/*... vtk elmentos*/
  #define VTK_TRIA      5
  #define VTK_QUAD      9
  #define VTK_TETR     10
  #define VTK_HEXA     12
/*...................................................................*/

/*... definicao do tipo de inteiros usados*/
  #define INT         int 
  #define INTC       "int"
  #define DOUBLE    double
  #define DOUBLEC  "double"
  #define LDOUBLE   long double
  #define LDOUBLEC  "long double"
  #ifdef _MSC_VER
  /*  #define LONG_INT __int64*/
    typedef __int64 LONG_INT;
    #define RESTRICT __restrict
  #else
    #define LONG_INT long  
  /*  typedef __int64 LONG_INT;*/
    #define RESTRICT restrict
  #endif
/*...................................................................*/

/*... macro para acesso matricial em vetores*/
  #define   MAT2D(i,j,vector,col)           vector[(i)*(col)+(j)]
  #define   MAT3D(i,j,k,vector,col1,col2)   vector[(i)*(col1)*(col2)+(col2)*(j)+(k)]
/*...................................................................*/

/*...*/
  #define   SKYLINEU(i,j,jDiag,vector)      vector[jDiag[(j)]+(i)-(j)]         
  #define   SKYLINEL(i,j,jDiag,vector)      vector[jDiag[(i)]+(j)-(i)]         
  #define   INDEXSKYLINEU(i,j,jDiag)        (jDiag[(j)]+(i)-(j))         
  #define   INDEXSKYLINEL(i,j,jDiag)        (jDiag[(i)]+(j)-(i))     
/*...................................................................*/

/*...*/
  #define HccaAbs(a)  (((a) < 0.0e0) ? (-a) : (a))
/*...................................................................*/

/*... definicao de funcoes*/
  #define min(a, b)  (((a) < (b)) ? (a) : (b))
  #define max(a, b)  (((a) > (b)) ? (a) : (b))
  #define vectorPlusOne(v,n,i)  for((i)=0;(i)<(n);(i++)) (v[(i)]++) 
  #define vectorMinusOne(v,n,i) for((i)=0;(i)<(n);(i++)) (v[(i)]--)  
/*...................................................................*/

/*...*/
  #define LDMT           1
  #define LDLT           2
  #define GGT            3
  #define LUKIJ          4
  #define LUIKJ          5
  #define LUKJI          6
  #define LUJKI          7
  #define DOOLITTLEL     8
  #define DOOLITTLELC    9
  #define DOOLITTLECL   10
  #define LUKIJPP       11
  #define DOOLITTLECLPP 12
/*...................................................................*/

/*... */            
  #define KELVINCONV 273.15E+00
  #define CELSIUS_FOR_KELVIN(t) ((t)+(KELVINCONV)) 
  #define KELVIN_FOR_CELSIUS(t) ((t)-(KELVINCONV)) 
  #define TEMP(tc,t,fKelvin)  if((fKelvin)){ (tc) = (t);}\
                              else{(tc) = CELSIUS_FOR_KELVIN((t));}
/*...................................................................*/

/*... Saida de Erro*/                                                  
  #define ERRO_RCM fprintf(stderr,"\nrcm - fatal error!\n")

  #define ERRO_OP(line,file,func,op)\
    fprintf(stderr,"Opecao %d e invalida!!\n",op);\
    fprintf(stderr,"Arquivo:%s\nFonte:  %s\nLinha:  %d\n"\
          ,file,func,line);\
    exit(EXIT_FAILURE);
  
  #define ERRO_NORM(line,file,func,lin,col)\
    fprintf(stderr,"Erro no calulo da norma!!\n");\
    fprintf(stderr,"Numero de linha  : %d\n"\
                   "Numero de colunas: %d\n",lin,col);\
    fprintf(stderr,"Arquivo:%s\nFonte:  %s\nLinha:  %d\n"\
          ,file,func,line);\
    exit(EXIT_FAILURE);
 
  #define ERRO_GERAL(file,func,str)\
    fprintf(stderr,"Erro: %s!!\n",str);\
    fprintf(stderr,"Arquivo:%s\nFonte:  %s\n",file,func);\
    exit(EXIT_FAILURE);

  #define ERRO_MALLOC(point,str,line,file,func)\
     if(point == NULL){\
     fprintf(stderr,"Erro na alocacao do vetor %s\n",str);\
     fprintf(stderr,"Arquivo:%s\nFonte:  %s\nLinha:  %d\n"\
            ,file,func,line);\
     exit(EXIT_FAILURE);}
/*...................................................................*/

/*...................................................................*/
#endif/*_DEFINE_*/
