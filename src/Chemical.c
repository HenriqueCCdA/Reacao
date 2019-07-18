#include<Chemical.h>

DOUBLE molarMass(DOUBLE *RESTRICT y     ,DOUBLE *RESTRICT mW
                ,short unsigned const nSp)
{
  short unsigned i;
  DOUBLE molarMass;

/*... massa especifica*/
  for(i=0,molarMass=0.e0;i<nSp;i++)
    molarMass += y[i]/mW[i];
  molarMass = 1.e0/molarMass;

  return molarMass;
}


/*********************************************************************
 * Data de criacao    : 19/06/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 * ------------------------------------------------------------------*
 * reverseK: calculo do coeficiente reverso atraves do direto        *
 * ------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 * ------------------------------------------------------------------* 
 * T     -> temperatura em kelvin                                    *
 * Rgas  -> constante dos gases ideias                               *
 * kf    -> coeficiente de reacao direto                             *
 * nSp   -> numero de especies                                       *
 * ------------------------------------------------------------------*
 * Paramanetros de saida:                                            *
 * ------------------------------------------------------------------*
 * ------------------------------------------------------------------*
 * OBS:                                                              *
 * ------------------------------------------------------------------*  
 *********************************************************************/
static DOUBLE reverseK(Reaction *reac ,PropVar *pGas 
                      ,DOUBLE const T ,DOUBLE const Rgas
                      ,DOUBLE const kf,DOUBLE const pres
                      ,unsigned short const nSp)
{
  unsigned short i;
  DOUBLE KP,de,sCoef,S,H,dG,tmp1,tmp2,tmp3=Rgas*T;

  for(i=0,sCoef=0.e0,dG=0.e0;i<nSp;i++)
  {
    de    = reac->stch[i][1] - reac->stch[i][0];
    sCoef+= de;
    H = polNasaH(pGas->pol+i,T,true);
    S = polNasaS(pGas->pol+i,T,true);
    dG+= de*(H - T*S);  
  }
  tmp1 = pow(pres/tmp3,sCoef);
  tmp2 = dG/tmp3;
  
  KP = tmp1*exp(-tmp2);
  return kf/KP;

}
/********************************************************************/

/*********************************************************************
 * Data de criacao    : 19/06/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 * ------------------------------------------------------------------*
 * arrheniusLaw: Lei de arrhenius                                    *
 * ------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 * ------------------------------------------------------------------* 
 * T     -> temperatura em kelvin                                    *
 * ------------------------------------------------------------------*
 * Paramanetros de saida:                                            *
 * ------------------------------------------------------------------*
 * ------------------------------------------------------------------*
 * OBS:                                                              *
 * ------------------------------------------------------------------*  
 *********************************************************************/
static DOUBLE arrheniusLaw(Arrhenius *arr,DOUBLE const T)
{
  DOUBLE A = arr->A, beta = arr->beta, Ta = arr->Ta;   

  return A*pow(T,beta)*exp(-Ta/T);    

}
/********************************************************************/

/*********************************************************************
 * Data de criacao    : 19/06/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 * ------------------------------------------------------------------*
 * massActionMass:                                                   *
 * ------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 * ------------------------------------------------------------------* 
 * c     -> concentracao molar das especies                          *
 * w     -> nao definido                                             *
 * T     -> temperatura em kelvin                                    *
 * nSp   -> numero de especies                                       *
 * ------------------------------------------------------------------*
 * Paramanetros de saida:                                            *
 * ------------------------------------------------------------------*
 * ------------------------------------------------------------------*
 * OBS:                                                              *
 * ------------------------------------------------------------------*  
 *********************************************************************/
static DOUBLE massActionMass(Reaction *reac,PropVar *pGas 
                     ,DOUBLE  *RESTRICT c
                     ,DOUBLE const T,unsigned short const nSp)
{
  unsigned short i;
  DOUBLE kf,kr,qf,qr;  

/*... reacao direta*/
  kf = arrheniusLaw(&reac->ArrF, T);

  for(i=0,qf=kf;i<nSp;i++)
    qf *= pow(c[i],reac->exp[i][0]);

/*... reacao inversao*/
  qr = 0.e0;
  if(reac->reverse)
  {
/*  kr = reverseK(reac       ,pGas
                 ,T          ,R
                 ,kf         ,1.0     
                 ,nSp);*/
    kr = arrheniusLaw(&reac->ArrR, T);

    for(i=0,qr=kr;i<nSp;i++)
      qr *= pow(c[i],reac->exp[i][1]);
  }

  return qf-qr;
}
/********************************************************************/

/*********************************************************************
 * Data de criacao    : 19/06/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 * ------------------------------------------------------------------*
 * molaRateReaction: taxa de reacao molar                            *
 * ------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 * ------------------------------------------------------------------* 
 * c     -> concentracao molar das especies                          *
 * w     -> nao definido                                             *
 * T     -> temperatura em kelvin                                    *
 * ------------------------------------------------------------------*
 * Paramanetros de saida:                                            *
 * ------------------------------------------------------------------*
 * w     -> taxa de reacao das especies                              *
 * ------------------------------------------------------------------*
 * OBS:                                                              *
 * ------------------------------------------------------------------*  
 *********************************************************************/
static void molarRateReaction(Chemical *chem,PropVar *pGas
                      ,DOUBLE *RESTRICT c,DOUBLE *RESTRICT w
                      ,DOUBLE const T    )
{
  unsigned short i, j, nReac =chem->nReac, nSp = chem->nSp;
  DOUBLE q,de;
  
  for(j=0;j<nSp;j++)
    w[j] = 0.e0;

  for(i=0;i<nReac;i++)
  {
    q = massActionMass(chem->reac+i,pGas,c,T,nSp);    
    for(j=0;j<nSp;j++)
    {
/*... */
      de    = chem->reac[i].stch[j][1] - chem->reac[i].stch[j][0];
/*... kmol/m3s*/
      w[j] += 1.0e+03*de*q;
    }
  }
}
/********************************************************************/

/*********************************************************************
 * Data de criacao    : 19/06/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 * ------------------------------------------------------------------*
 * energySource: taxa de liberação de energia                        *
 * ------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 * ------------------------------------------------------------------* 
 * pVar  -> propriedades variaveis                                   *
 * w     -> taxa de reacao das especies                              *
 * mW    -> massa molar das especies                                 *
 * T     -> temperatura em kelvin                                    *
 * nSp   -> numero das especies                                      *
 * ------------------------------------------------------------------*
 * Paramanetros de saida:                                            *
 * ------------------------------------------------------------------*
 * ------------------------------------------------------------------*
 * OBS:                                                              *
 * ------------------------------------------------------------------*  
 *********************************************************************/
static DOUBLE energySource(PropVar *pVar
                   ,DOUBLE *RESTRICT w,DOUBLE *RESTRICT mW
                   ,DOUBLE const T    ,unsigned short nSp)
{
  unsigned short i;
  DOUBLE h,wT;  

  for(i=0,wT=0.e0;i<nSp;i++)
  {
    h = polNasaH(pVar->pol+i,T,true);
    wT += h*w[i];
  }
  
  return -wT;
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 19/06/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 * ------------------------------------------------------------------*
 * yFracSum: soma das fracoes massicas                               *
 * ------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 * ------------------------------------------------------------------*  
 * yFrac -> fracoes massicas                                         *
 * n     -> numero de especies                                       * 
 * ------------------------------------------------------------------*
 * Paramanetros de saida:                                            *
 * ------------------------------------------------------------------*
 * ------------------------------------------------------------------*
 * OBS:                                                              *
 * ------------------------------------------------------------------*  
 * *******************************************************************/
DOUBLE yFracSum(DOUBLE *RESTRICT yFrac, unsigned short n) 
{
  short i;
  DOUBLE s;
  
  for(i=0,s=0.e0;i<n;i++)
    s += yFrac[i];

  return s;
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 19/06/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 * ------------------------------------------------------------------*
 * chemicalRsh :                                                     *
 * ------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 * ------------------------------------------------------------------* 
 * t  -> tempo das simulacao                                         *
 * y  -> solução                                                     *
 * w  -> nao definido                                                *
 * pt -> paramentro de calculo                                       *     
 * ------------------------------------------------------------------*
 * Paramanetros de saida:                                            *
 * ------------------------------------------------------------------*
 * w  -> vetor independente                                          *
 * ------------------------------------------------------------------*
 * OBS:                                                              *
 * ------------------------------------------------------------------*
 * y[0,1,2,...N-1] -> fracoes massicas                               *
 * y[N]            -> temperatura                                    *    
 * *******************************************************************/
void chemicalRsh(DOUBLE const t    ,DOUBLE *RESTRICT y
                ,DOUBLE *RESTRICT w,void **pt)
{
  unsigned short i,nSp;
  DOUBLE c[MAX_SPECIES],T;
  DOUBLE cp,cpk,rho,molarMassGas; 
  Chemical *chem = pt[0];
  PropVar  *pVar = pt[1];
  
  nSp = chem->nSp;
  T = y[nSp];

/*... massa especifica*/
  molarMassGas = molarMass(y,chem->mW,nSp);
  rho = (molarMassGas*PREREF)/(1.e+03*R*T);

/*... calor especifico*/
  for(i=0,cp=0.e0;i<nSp;i++)
  {
    cpk = polNasaCp(pVar->pol+i,T);
    cp += max(y[i],0.e0)*cpk;
  }

/*... converter para mol/cm^3*/
  for(i=0;i<nSp;i++)
    c[i] = 1.e-03*rho*max(y[i],0.e0)/chem->mW[i];
/*...................................................................*/

  molarRateReaction(chem,pVar,c,w,T);
  w[nSp] = energySource(pVar,w,chem->mW,T,nSp);
  
  for(i=0;i<nSp;i++)
    w[i] *= (chem->mW[i]/rho);
  w[nSp] /= (cp*rho);
}
/*********************************************************************/