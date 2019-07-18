#include"Edo.h"

/*********************************************************************
 * Data de criacao    : 16/06/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 * ------------------------------------------------------------------*
 * stepperRossControllerSuccess :                                    *
 * ------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 * ------------------------------------------------------------------*
 * ------------------------------------------------------------------*
 * Paramanetros de saida:                                            *
 * ------------------------------------------------------------------*
 * ------------------------------------------------------------------*
 * OBS:                                                              *
 * *******************************************************************/
static bool stepperRossControllerSuccess(DOUBLE const err, DOUBLE *h
                                        ,DOUBLE *errOld  , DOUBLE *hOld
                                        ,DOUBLE *hNext   , bool *reject
                                        ,bool *firstStep)
{
  static const DOUBLE safe = 0.9
                    , fac1 = 5.0
                    , fac2 = 1.0 / 6.0;
  DOUBLE fac;
  DOUBLE hNew;
  DOUBLE facPred;

/*...*/
  fac  = max(fac2, min(fac1, pow(err, 0.25) / safe));
  hNew = *h / fac;
/*...................................................................*/

/*... erro aceitavel*/
  if (err <= 1.0)
  {
    if (!(*firstStep))
    {
      facPred = (*hOld / *h)*pow(err*err / *errOld, 0.25) / safe;
      facPred = min(fac2, min(fac1, facPred));
      fac = min(fac, facPred);
      hNew = *h / fac;
    }
    *firstStep = 0;
    *hOld      = *h;
    *errOld    = max(0.01, err);
    if (*reject)
      hNew = (*h >= 0.0 ? min(hNew, *h) : min(hNew, *h));
    *hNext = hNew;
    *reject = 0;
    return 1;
  }
/*... err demasido indo para nova tentativa com um novo h*/
  else
  {
    *h = hNew;
    *reject = 1;
    return 0;
  }
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 16/06/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 * ------------------------------------------------------------------*
 * stepperDopr5 :  Metodo Runge-kutta explicito de 5 ° ordem como    *
 * passo adaptativo                                                  *
 * ------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 * ------------------------------------------------------------------*
 * y     -> valores inicial                                          *
 * param -> paramatros utilizado no calculo da funcao f              *
 * x1    -> limite inferio do intervalo de integracao                *
 * x2    -> limite superior do intervalo de integraçao               *
 * hInit -> passo inicial                                            *
 * maxIt -> numero maximo de iteração no processo de integracao      *
 * nEdo  -> numero de edos                                           *
 * maxInt-> numero maximo de passos de integracao                    * 
 * fStopIt -> true para quando maxInt é atingido                     *
 *           false nao a limite de passo de initegra                 *
 * outPut  -> cod de escrita                                         *
 *            1 - tela                                               *
 *            2 - arquivo                                            *
 * fName   -> nomde do aquivo de saida                               *
 * rhs     -> f                                                      *
 * jacY    -> parte do jacabiano da f dependo explicidamente de y    *
 * jacX    -> parte do jacabiano da f dependo explicidamente de x    *
 * ------------------------------------------------------------------*
 * Paramanetros de saida:                                            *
 * ------------------------------------------------------------------*
 * ------------------------------------------------------------------*
 * *******************************************************************/
int stepperRoss(DOUBLE *y                 , DOUBLE *param
               , DOUBLE const x1          , DOUBLE const x2
               , DOUBLE const hInit       , DOUBLE const aTol
               , DOUBLE const rTol        , short const nEdo
               , unsigned INT const maxInt, bool const fStopIt
               , short const  outCod      , const char *const fName
               , void (*rhs)()
               , void (*jacY)()           , void (*jacX)())
{
/*... constantes*/
  DOUBLE const c2 = 0.386;
  DOUBLE const c3 = 0.21;
  DOUBLE const c4 = 0.63;
  DOUBLE const bet2p = 0.0317;
  DOUBLE const bet3p = 0.0635;
  DOUBLE const bet4p = 0.3438;
  DOUBLE const d1 = 0.2500000000000000e+00;
  DOUBLE const d2 = -0.1043000000000000e+00;
  DOUBLE const d3 = 0.1035000000000000e+00;
  DOUBLE const d4 = -0.3620000000000023e-01;
  DOUBLE const a21 = 0.1544000000000000e+01;
  DOUBLE const a31 = 0.9466785280815826e+00;
  DOUBLE const a32 = 0.2557011698983284e+00;
  DOUBLE const a41 = 0.3314825187068521e+01;
  DOUBLE const a42 = 0.2896124015972201e+01;
  DOUBLE const a43 = 0.9986419139977817e+00;
  DOUBLE const a51 = 0.1221224509226641e+01;
  DOUBLE const a52 = 0.6019134481288629e+01;
  DOUBLE const a53 = 0.1253708332932087e+02;
  DOUBLE const a54 = -0.6878860361058950e+00;
  DOUBLE const c21 = -0.5668800000000000e+01;
  DOUBLE const c31 = -0.2430093356833875e+01;
  DOUBLE const c32 = -0.2063599157091915e+00;
  DOUBLE const c41 = -0.1073529058151375e+00;
  DOUBLE const c42 = -0.9594562251023355e+01;
  DOUBLE const c43 = -0.2047028614809616e+02;
  DOUBLE const c51 = 0.7496443313967647e+01;
  DOUBLE const c52 = -0.1024680431464352e+02;
  DOUBLE const c53 = -0.3399990352819905e+02;
  DOUBLE const c54 = 0.1170890893206160e+02;
  DOUBLE const c61 = 0.8083246795921522e+01;
  DOUBLE const c62 = -0.7981132988064893e+01;
  DOUBLE const c63 = -0.3152159432874371e+02;
  DOUBLE const c64 = 0.1631930543123136e+02;
  DOUBLE const c65 = -0.6058818238834054e+01;
  DOUBLE const gam = 0.2500000000000000e+00;
  DOUBLE const d21 = 0.1012623508344586e+02;
  DOUBLE const d22 = -0.7487995877610167e+01;
  DOUBLE const d23 = -0.3480091861555747e+02;
  DOUBLE const d24 = -0.7992771707568823e+01;
  DOUBLE const d25 = 0.1025137723295662e+01;
  DOUBLE const d31 = -0.6762803392801253e+00;
  DOUBLE const d32 = 0.6087714651680015e+01;
  DOUBLE const d33 = 0.1643084320892478e+02;
  DOUBLE const d34 = 0.2476722511418386e+02;
  DOUBLE const d35 = -0.6594389125716872e+01;
/*.....................................................................*/
  FILE *file = NULL;
  bool accept, reject, firstStep;
  short i, j, maxChangeSteps, iChangeSteps;
  unsigned INT it,IntegralStepMax;
  INT *p = NULL;
  DOUBLE *a     = NULL
        ,*f     = NULL
        ,*fNew  = NULL
        ,*dfdy  = NULL
        ,*dfdx  = NULL
        ,*yt    = NULL
        ,*yOut  = NULL
        ,*yErr  = NULL
        ,*g1    = NULL
        ,*g2    = NULL
        ,*g3    = NULL
        ,*g4    = NULL
        ,*g5    = NULL
        ,*g6    = NULL
        ,*w     = NULL;

  DOUBLE x,h,xph,err,errOld, hOld, hNext;

/*...*/
  maxChangeSteps  = 100;
  IntegralStepMax = maxInt;
/*...................................................................*/

/*...*/
  if (outCod == FILE_OUT)
    file = openFile(fName, "w");
/*...................................................................*/

/*...*/
  a    = (DOUBLE*)malloc(nEdo * nEdo * sizeof(DOUBLE));
  testAlloc(a,"a");
  dfdy = (DOUBLE*)malloc(nEdo * nEdo * sizeof(DOUBLE));
  testAlloc(dfdy,"dfdy");
  dfdx = (DOUBLE*)malloc(nEdo * sizeof(DOUBLE));
  testAlloc(dfdx,"dfdx");
  f    = (DOUBLE*)malloc(nEdo * sizeof(DOUBLE));
  testAlloc(f,"f");
  fNew = (DOUBLE*)malloc(nEdo * sizeof(DOUBLE));
  testAlloc(fNew,"fNew");
  yt   = (DOUBLE*)malloc(nEdo * sizeof(DOUBLE));
  testAlloc(yt,"yt");
  yOut = (DOUBLE*)malloc(nEdo * sizeof(DOUBLE));
  testAlloc(yOut,"yOut");
  yErr = (DOUBLE*)malloc(nEdo * sizeof(DOUBLE));
  testAlloc(yErr,"yErr");
  g1   = (DOUBLE*)malloc(nEdo * sizeof(DOUBLE));
  testAlloc(g1,"g1");
  g2   = (DOUBLE*)malloc(nEdo * sizeof(DOUBLE));
  testAlloc(g2,"g2");
  g3   = (DOUBLE*)malloc(nEdo * sizeof(DOUBLE));
  testAlloc(g3,"g3");
  g4   = (DOUBLE*)malloc(nEdo * sizeof(DOUBLE));
  testAlloc(g4,"g4");
  g5   = (DOUBLE*)malloc(nEdo * sizeof(DOUBLE));
  testAlloc(g5,"g5");
  g6   = (DOUBLE*)malloc(nEdo * sizeof(DOUBLE));
  testAlloc(g6,"g6");
  w    = (DOUBLE*)malloc(nEdo * sizeof(DOUBLE));
  testAlloc(w,"w");
  p    = (INT   *)malloc(nEdo * sizeof(INT));
  testAlloc(p,"p");
/*...................................................................*/

/*...*/
  h            = hInit;
  it           = 0;
  x            = x1;
  errOld       = 1.e0;
/*...................................................................*/

/*...*/
  firstStep = true;
  do
  {
/*... tentativa com o passo de tamanho h*/
    reject       = false;
    iChangeSteps = false;
/*... derivada df/dy*/
    jacY(x,y, dfdy,param);
/*... derivada df/dx*/
    jacX(x, y, dfdx, param);
/*... f(y)*/
    rhs(x, y, f, param);
/*...*/
    do 
    {
/*... A = (1/gammah - dydx)*/
      for (i = 0; i<nEdo; i++)
      { 
        for (j = 0; j<nEdo; j++) 
          MAT2D(i,j,a,nEdo) = -MAT2D(i, j, dfdy, nEdo);
        MAT2D(i, i, a, nEdo) += 1.e0 / (gam*h);
      }
/*...................................................................*/

/*... decomposicao de a*/
      fatLUpp(a,p, nEdo, LUKIJPP);
/*...................................................................*/
  
/*...*/
      for (i = 0; i<nEdo; i++)
        yt[i] = f[i] + h * d1*dfdx[i];
/*...................................................................*/

/*... solver g1*/
      solvLUpp(a, yt, w, p, nEdo);
      for (i = 0; i<nEdo; i++)
      {
        g1[i] = yt[i];
        yt[i] = y[i] + a21 * g1[i];
      }
/*...................................................................*/

/*...*/
      rhs(x + c2 * h, yt, fNew, param);
      for (i = 0; i<nEdo; i++)
        yt[i] = fNew[i] + h * d2*dfdx[i] + c21 * g1[i] / h;
/*...................................................................*/

/*... solver g2*/
      solvLUpp(a, yt, w, p, nEdo);
      for (i = 0; i<nEdo; i++)
      {
        g2[i] = yt[i];
        yt[i] = y[i] + a31 * g1[i] + a32 * g2[i];
      }
/*...................................................................*/

/*...*/
      rhs(x + c3 * h, yt, fNew,param);
      for (i = 0; i<nEdo; i++)
        yt[i] = fNew[i] + h * d3*dfdx[i] + (c31*g1[i] + c32*g2[i]) / h;
/*...................................................................*/

/*... solver g3*/
      solvLUpp(a, yt, w, p, nEdo);
      for (i = 0; i<nEdo; i++)
      {
        g3[i] = yt[i];
        yt[i] = y[i] + a41*g1[i] + a42*g2[i] + a43*g3[i];
      }
/*...................................................................*/

/*...*/
      rhs(x + c4 * h, yt, fNew, param);
      for (i = 0; i<nEdo; i++)
        yt[i] = fNew[i] + h * d4*dfdx[i] + (c41*g1[i] + c42*g2[i] 
                                                    + c43*g3[i]) / h;
/*...................................................................*/

/*... solver g4*/
      solvLUpp(a, yt, w, p, nEdo);
      for (i = 0; i<nEdo; i++)
      {
        g4[i] = yt[i];
        yt[i] = y[i] + a51*g1[i] + a52*g2[i] + a53*g3[i] + a54*g4[i];
      }
/*...................................................................*/

/*...*/
      xph = x + h;
      rhs(xph, yt, fNew, param);
      for (i = 0; i<nEdo; i++)
        g6[i] = fNew[i] + (c51*g1[i] + c52*g2[i] + c53*g3[i] 
                                                 + c54*g4[i]) / h;
/*...................................................................*/

/*... solver g5*/
      solvLUpp(a, g6, w, p, nEdo);
      for (i = 0; i<nEdo; i++)
      {
        g5[i]  = g6[i];
        yt[i] += g5[i];
      }
/*...................................................................*/

/*...*/
      rhs(xph, yt, fNew, param);
      for (i = 0; i<nEdo; i++)
        g6[i] = fNew[i] + (c61*g1[i] + c62*g2[i] + c63*g3[i] 
                         + c64*g4[i] + c65*g5[i]) / h;
/*...................................................................*/

/*... solver g6*/
      solvLUpp(a, g6, w, p, nEdo);
      for (i = 0; i<nEdo; i++)
      {
        yErr[i] = g6[i];
        yOut[i] = yt[i] + yErr[i];
      }
/*...................................................................*/

/*... calculo do erro*/
      err = edoError(y   , yOut, yErr
                     ,aTol, rTol, nEdo);
/*...................................................................*/

/*... teste a solucao para o h*/
      accept = stepperRossControllerSuccess(err         , &h
                                     , &errOld    , &hOld
                                     , &hNext     , &reject
                                     , &firstStep);

      iChangeSteps++;
    }while(!accept && iChangeSteps < maxChangeSteps);
/*...................................................................*/

/*...*/
    if(iChangeSteps == maxChangeSteps)
    {
      printf("Numero maximo de mudanca de passo alcancada!!\n");
      exit(EXIT_FAILURE);
    }
/*...................................................................*/

/*...*/
    for (i = 0; i<nEdo; i++)
      y[i] = yOut[i];
/*...................................................................*/

/*...*/
     writeOutput(it            , x
                ,h             , y
                ,nEdo          , outCod
                ,file);
/*.....................................................................*/

/*...*/
    it++;
    if (x == x + h)
    {
      printf("stepsize underflow in stepperRoss\n");
      exit(EXIT_FAILURE);
    }
    x += h;
    h = hNext;
/*.....................................................................*/

/*...*/
    if(it > IntegralStepMax && fStopIt)
    {
       printf("Numero maximo de passos para integracao alcancado !!\n");
      exit(EXIT_FAILURE);
    }
/*.....................................................................*/
  } while ((x - h) <= x2);
/*.....................................................................*/

/*...*/
  free(a   );
  free(dfdy);
  free(dfdx);
  free(f   );
  free(fNew);
  free(yt  );
  free(yOut);
  free(yErr);
  free(g1  );
  free(g2  );
  free(g3  );
  free(g4  );
  free(g5  );
  free(g6  );
/*.....................................................................*/

  if(outCod == FILE_OUT)
    fclose(file);

  return it;
}
/***********************************************************************/