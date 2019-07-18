#include<string.h>

#include<Edo.h>
#include<NasaPol.h>
#include<StructDefine.h>
#include<Chemical.h>

void setPolNasa(PropVar *pGas,DOUBLE *mW,short const nSp);

void stchConvMass(Chemical *chem);

int main(int argc, char *argv[])
{
  unsigned short i;
  Chemical chem;
  PropVar pGas;
  INT nStep;
  DOUBLE y[MAX_SPECIES+1];
  DOUBLE tol,h,t0,tF;
  void *pt[10];
/*...*/
  pt[0]      =&chem;
  pt[1]      =&pGas;  
  chem.nSp   = 6;
  chem.nReac = 2;
/*.....................................................................*/

  strcpy(chem.name[0],"CH4");
  chem.mW[0] = MW_C +  4*MW_H;

  strcpy(chem.name[1],"O2");
  chem.mW[1] = 2*MW_O;

  strcpy(chem.name[2],"H2O");
  chem.mW[2] = 2*MW_H +  MW_O;

  strcpy(chem.name[3],"CO2");
  chem.mW[3] = MW_C +  2*MW_O;

  strcpy(chem.name[4],"CO");
  chem.mW[4] = MW_C +  MW_O;

  strcpy(chem.name[5],"N2");
  chem.mW[5] = 2*MW_N;
/*... CH4 + O2 -> CO + 2H2O*/

/*... reagente*/
  chem.reac[0].stch[0][0] = 1.0;
  chem.reac[0].stch[1][0] = 1.5;
  chem.reac[0].stch[2][0] = 0.0;
  chem.reac[0].stch[3][0] = 0.0;
  chem.reac[0].stch[4][0] = 0.0;
  chem.reac[0].stch[5][0] = 0.0;
/*... produto*/
  chem.reac[0].stch[0][1] = 0.0;
  chem.reac[0].stch[1][1] = 0.0;
  chem.reac[0].stch[2][1] = 2.0;
  chem.reac[0].stch[3][1] = 0.0;
  chem.reac[0].stch[4][1] = 1.0;
  chem.reac[0].stch[5][1] = 0.0;
  

  chem.reac[0].reverse     = false; 
  chem.reac[0].ArrF.A      = 2.00E+15;
  chem.reac[0].ArrF.beta   = 0.0;
  chem.reac[0].ArrF.E      = 35000.0;
  chem.reac[0].ArrF.Ta     = chem.reac[0].ArrF.E/RC;
  chem.reac[0].exp[0][0]   = 0.9e0;
  chem.reac[0].exp[1][0]   = 1.1e0;
  chem.reac[0].exp[2][0]   = 0.0e0;
  chem.reac[0].exp[3][0]   = 0.0e0;
  chem.reac[0].exp[4][0]   = 0.0e0;
  chem.reac[0].exp[5][0]   = 0.0e0;

  chem.reac[0].exp[0][1]   = 0.0e0;
  chem.reac[0].exp[1][1]   = 0.0e0;
  chem.reac[0].exp[2][1]   = 0.0e0;
  chem.reac[0].exp[3][1]   = 0.0e0;
  chem.reac[0].exp[4][1]   = 0.0e0;
  chem.reac[0].exp[5][1]   = 0.0e0;
/*.....................................................................*/

/*... CO + 1/2 O2 -> CO2*/
/*... reagente*/
  chem.reac[1].stch[0][0] = 0.0;
  chem.reac[1].stch[1][0] = 0.5;
  chem.reac[1].stch[2][0] = 0.0;
  chem.reac[1].stch[3][0] = 0.0;
  chem.reac[1].stch[4][0] = 1.0;
  chem.reac[1].stch[5][0] = 0.0;
/*... produto*/
  chem.reac[1].stch[0][1] = 0.0;
  chem.reac[1].stch[1][1] = 0.0;
  chem.reac[1].stch[2][1] = 0.0;
  chem.reac[1].stch[3][1] = 1.0;
  chem.reac[1].stch[4][1] = 0.0;
  chem.reac[1].stch[5][1] = 0.0;

  chem.reac[1].reverse    = true; 
  chem.reac[1].ArrF.A     = 2.00E+09;
  chem.reac[1].ArrF.beta  = 0.0;
  chem.reac[1].ArrF.E     = 12000.0;
  chem.reac[1].ArrF.Ta    = chem.reac[1].ArrF.E/RC;

  chem.reac[1].ArrR.A     = 8.11704E+10;
  chem.reac[1].ArrR.beta  = 0.0;
  chem.reac[1].ArrR.E     = 77194.0;
  chem.reac[1].ArrR.Ta    = chem.reac[1].ArrR.E/RC;

  chem.reac[1].exp[0][0]  = 0.00;
  chem.reac[1].exp[1][0]  = 0.50;
  chem.reac[1].exp[2][0]  = 0.00;
  chem.reac[1].exp[3][0]  = 0.00;
  chem.reac[1].exp[4][0]  = 1.00;
  chem.reac[1].exp[5][0]  = 0.00;

  chem.reac[1].exp[0][1]  = 0.00;
  chem.reac[1].exp[1][1]  = 0.00;
  chem.reac[1].exp[2][1]  = 0.00;
  chem.reac[1].exp[3][1]  = 1.00;
  chem.reac[1].exp[4][1]  = 0.00;
  chem.reac[1].exp[5][1]  = 0.00;
/*.....................................................................*/

/*... testa os coeficiente estequiometricos*/
   stchConvMass(&chem);
/*.....................................................................*/

/*...*/
  t0          =  0.0e0;
  tF          =  0.01e0;
  y[chem.nSp] = 800e0;
  y[0]        = 0.1e0;
  y[1]        = 0.6e0;
  y[2]        = 0.e0;
  y[3]        = 0.e0;
  y[4]        = 0.e0;
  y[5]        = 0.3e0;
/*.....................................................................*/

/*...*/
  setPolNasa(&pGas,chem.mW,chem.nSp);
/*.....................................................................*/

/*...*/
  h   = 1.e-06;
  tol = 1.e-11;
/*.....................................................................*/

/*...*/
  nStep = StepperSie(y          , pt   
           , t0                 , tF
           , h                  , 1.0e-5
           , tol                , tol   
           , chem.nSp + 1
           , 2000000            , true
           , FILE_OUT           , "teste.out"
           , chemicalRsh        , NULL); 
/*.....................................................................*/

  printf("nStep : %d\n",nStep);
  for(i=0;i<chem.nSp;i++)
    printf("%-20s : %.4lf\n",chem.name[i],y[i]);
  printf("%-20s : %lf\n","T(K)",y[chem.nSp]);
  printf("%-20s : %lf\n","sY",yFracSum(y,chem.nSp));

  return EXIT_SUCCESS;

}

void setPolNasa(PropVar *pGas,DOUBLE *mW,short const nSp)
{
  unsigned short i,j,k;
  DOUBLE h;
/*... pol nasa*/

/*...CH4*/
  pGas->pol[0].type        = 7;
  pGas->pol[0].range[0][0] = 200.e0;
  pGas->pol[0].range[0][1] = 1000.e0;
  pGas->pol[0].a1[0][0]    = 5.14825732E+00;
  pGas->pol[0].a1[0][1]    =-1.37002410E-02;
  pGas->pol[0].a1[0][2]    = 4.93749414E-05;
  pGas->pol[0].a1[0][3]    =-4.91952339E-08;
  pGas->pol[0].a1[0][4]    = 1.70097299E-11;
  pGas->pol[0].a1[0][5]    =-1.02453222E+04;
  pGas->pol[0].a1[0][6]    =-4.63322726E+00; 
             
  pGas->pol[0].range[1][0] = 1000.e0;
  pGas->pol[0].range[1][1] = 6000.e0;
  pGas->pol[0].a1[1][0]    = 1.91178600E+00;
  pGas->pol[0].a1[1][1]    = 9.60267960E-03;
  pGas->pol[0].a1[1][2]    =-3.38387841E-06;
  pGas->pol[0].a1[1][3]    = 5.38797240E-10;
  pGas->pol[0].a1[1][4]    =-3.19306807E-14;
  pGas->pol[0].a1[1][5]    =-1.00992136E+04;
  pGas->pol[0].a1[1][6]    = 8.48241861E+00; 
 
/*...O2*/
  pGas->pol[1].type        = 7;
  pGas->pol[1].range[0][0] = 200.e0;
  pGas->pol[1].range[0][1] = 1000.e0;
  pGas->pol[1].a1[0][0]    = 3.78245636E+00;
  pGas->pol[1].a1[0][1]    =-2.99673415E-03;
  pGas->pol[1].a1[0][2]    = 9.84730200E-06;
  pGas->pol[1].a1[0][3]    =-9.68129508E-09;
  pGas->pol[1].a1[0][4]    = 3.24372836E-12;
  pGas->pol[1].a1[0][5]    =-1.06394356E+03;
  pGas->pol[1].a1[0][6]    = 3.65767573E+00;
             
  pGas->pol[1].range[1][0] = 1000.e0;
  pGas->pol[1].range[1][1] = 6000.e0;
  pGas->pol[1].a1[1][0]    = 3.66096083E+00;
  pGas->pol[1].a1[1][1]    = 6.56365523E-04;
  pGas->pol[1].a1[1][2]    =-1.41149485E-07;
  pGas->pol[1].a1[1][3]    = 2.05797658E-11;
  pGas->pol[1].a1[1][4]    =-1.29913248E-15;
  pGas->pol[1].a1[1][5]    =-1.21597725E+03;
  pGas->pol[1].a1[1][6]    = 3.41536184E+00; 

/*... H2O*/
  pGas->pol[2].type        = 7;
  pGas->pol[2].range[0][0] = 200.e0;
  pGas->pol[2].range[0][1] = 1000.e0;
  pGas->pol[2].a1[0][0]    = 4.19863520E+00;
  pGas->pol[2].a1[0][1]    =-2.03640170E-03;
  pGas->pol[2].a1[0][2]    = 6.52034160E-06;
  pGas->pol[2].a1[0][3]    =-5.48792690E-09;
  pGas->pol[2].a1[0][4]    = 1.77196800E-12;
  pGas->pol[2].a1[0][5]    =-3.02937260E+04;
  pGas->pol[2].a1[0][6]    =-8.49009010E-01;
             
  pGas->pol[2].range[1][0] = 1000.e0;
  pGas->pol[2].range[1][1] = 6000.e0;
  pGas->pol[2].a1[1][0]    = 2.67703890E+00;
  pGas->pol[2].a1[1][1]    = 2.97318160E-03;
  pGas->pol[2].a1[1][2]    =-7.73768890E-07;
  pGas->pol[2].a1[1][3]    = 9.44335140E-11; 
  pGas->pol[2].a1[1][4]    =-4.26899910E-15;
  pGas->pol[2].a1[1][5]    =-2.98858940E+04;
  pGas->pol[2].a1[1][6]    = 6.88255000E+00; 

/*... CO2*/
  pGas->pol[3].type        = 7;
  pGas->pol[3].range[0][0] = 200.e0;
  pGas->pol[3].range[0][1] = 1000.e0;
  pGas->pol[3].a1[0][0]    = 0.23568130E+01;
  pGas->pol[3].a1[0][1]    = 0.89841299E-02;
  pGas->pol[3].a1[0][2]    =-0.71220632E-05;
  pGas->pol[3].a1[0][3]    = 0.24573008E-08 ;
  pGas->pol[3].a1[0][4]    =-0.14288548E-12;
  pGas->pol[3].a1[0][5]    =-0.48371971E+05;
  pGas->pol[3].a1[0][6]    = 0.99009035E+01;
             
  pGas->pol[3].range[1][0] = 1000.e0;
  pGas->pol[3].range[1][1] = 6000.e0;
  pGas->pol[3].a1[1][0]    = 0.46365111E+01;
  pGas->pol[3].a1[1][1]    = 0.27414569E-02;
  pGas->pol[3].a1[1][2]    =-0.99589759E-06;
  pGas->pol[3].a1[1][3]    = 0.16038666E-09;
  pGas->pol[3].a1[1][4]    =-0.91619857E-14;
  pGas->pol[3].a1[1][5]    =-0.49024904E+05;
  pGas->pol[3].a1[1][6]    =-0.19348955E+01; 

/*...CO*/
  pGas->pol[4].type        = 7;
  pGas->pol[4].range[0][0] = 200.e0;
  pGas->pol[4].range[0][1] = 1000.e0;
  pGas->pol[4].a1[0][0]    = 3.57953350E+00;
  pGas->pol[4].a1[0][1]    =-6.10353690E-04;
  pGas->pol[4].a1[0][2]    = 1.01681430E-06;
  pGas->pol[4].a1[0][3]    = 9.07005860E-10;
  pGas->pol[4].a1[0][4]    =-9.04424490E-13;
  pGas->pol[4].a1[0][5]    =-1.43440860E+04;
  pGas->pol[4].a1[0][6]    = 3.50840930E+00;
             
  pGas->pol[4].range[1][0] = 1000.e0;
  pGas->pol[4].range[1][1] = 6000.e0;
  pGas->pol[4].a1[1][0]    = 3.04848590E+00;
  pGas->pol[4].a1[1][1]    = 1.35172810E-03;
  pGas->pol[4].a1[1][2]    =-4.85794050E-07;
  pGas->pol[4].a1[1][3]    = 7.88536440E-11;
  pGas->pol[4].a1[1][4]    =-4.69807460E-15;
  pGas->pol[4].a1[1][5]    =-1.42661170E+04;
  pGas->pol[4].a1[1][6]    = 6.01709770E+00; 

/*...N2*/
  pGas->pol[5].type        = 7;
  pGas->pol[5].range[0][0] = 200.e0;
  pGas->pol[5].range[0][1] = 1000.e0;
  pGas->pol[5].a1[0][0]    = 3.53100528E+00;
  pGas->pol[5].a1[0][1]    =-1.23660988E-04;
  pGas->pol[5].a1[0][2]    =-5.02999433E-07;
  pGas->pol[5].a1[0][3]    = 2.43530612E-09;
  pGas->pol[5].a1[0][4]    =-1.40881235E-12;
  pGas->pol[5].a1[0][5]    =-1.04697628E+03;
  pGas->pol[5].a1[0][6]    = 2.96747038E+00;
             
  pGas->pol[5].range[1][0] = 1000.e0;
  pGas->pol[5].range[1][1] = 6000.e0;
  pGas->pol[5].a1[1][0]    = 2.95257637E+00;
  pGas->pol[5].a1[1][1]    = 1.39690040E-03;
  pGas->pol[5].a1[1][2]    =-4.92631603E-07;
  pGas->pol[5].a1[1][3]    = 7.86010195E-11;
  pGas->pol[5].a1[1][4]    =-4.60755204E-15;
  pGas->pol[5].a1[1][5]    =-9.23948688E+02;
  pGas->pol[5].a1[1][6]    = 5.87188762E+00; 

/*...*/
  for(i=0;i<nSp;i++)
  {
    h = R/mW[i];
    for(k=0;k<2;k++)
      for(j=0;j<7;j++)
        pGas->pol[i].a2[k][j] = pGas->pol[i].a1[k][j]*h;
  }
/*.....................................................................*/

}

void stchConvMass(Chemical *chem)
{
  unsigned short i, j, nReac =chem->nReac, nSp = chem->nSp;
  DOUBLE de,mass;

/*...*/
  for(i=0;i<nReac;i++)
  {

    for(j=0,mass=0.e0;j<nSp;j++)
    {
/*... */
      de    = chem->reac[i].stch[j][1] - chem->reac[i].stch[j][0];
      mass += de*chem->mW[j];
    }
/*...................................................................*/

/*...*/
    if ( fabs(mass) > 1.e-15)
    {
      printf("Conservao de massa violada na equação quimica %hd %lf\n"
            ,i,mass);
    }
  }
/*...................................................................*/
}