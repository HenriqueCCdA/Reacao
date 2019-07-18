#ifndef _NASA_POL_H_
  #define _NASA_POL_H_

/*...*/
  #include<math.h>
/*...*/ 
  #include<Define.h> 
  #include<HccaStdBool.h>
/*...................................................................*/

/*...*/
   typedef struct
   {
     short type;
     DOUBLE range[3][2];
     DOUBLE a1[2][9],a2[2][9];  /*(KJ/KmolK) (KJ/KGK)*/
   }PolNasa;
/*...................................................................*/

  void nasaPolRange(PolNasa *a      , DOUBLE const x
                 ,DOUBLE **c      , DOUBLE *xNew
                 ,short const iCod);

  DOUBLE polNasaCp(PolNasa *a     , DOUBLE const x);

  DOUBLE polNasaH(PolNasa *a     , DOUBLE const x, bool const fKelvin);
  DOUBLE polNasaS(PolNasa *a     , DOUBLE const x, bool const fKelvin);

#endif /*_NASA_POL_H_*/