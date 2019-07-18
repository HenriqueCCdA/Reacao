#ifndef _EDO_H_
  #define _EDO_H_
/*...*/
  #include<stdlib.h>
  #include<stdio.h>
  #include<math.h>
/*...................................................................*/

/*...*/ 
  #include<Define.h> 
  #include<Gauss.h>
/*...................................................................*/

  #define SCREEN_OUT 1
  #define FILE_OUT   2
  
  FILE* openFile(const char* const name, const char* const mod);
  
  void rtk4(DOUBLE *y            , DOUBLE *param
           ,DOUBLE const x1      , DOUBLE const x2
           ,DOUBLE const hInit
           ,INT const maxIt     , short const nEdo
           ,short const outCod  , void (*rhs)());
  
  int stepperDopr5(DOUBLE *y                  , DOUBLE *param
                   ,DOUBLE const x1           , DOUBLE const x2
                   ,DOUBLE const aTol         , DOUBLE const rTol
                   ,DOUBLE const hInit        , short const nEdo
                   ,unsigned INT const maxInt , bool const fStopIt
                   ,const short outCod        , const char *const fName
                   ,void (*rhs)());
  
  int stepperDopr853(DOUBLE *y              , DOUBLE *param
                 ,DOUBLE const x1           , DOUBLE const x2
                 ,DOUBLE const hInit         ,DOUBLE const aTol         
                 , DOUBLE const rTol        , const short nEdo
                 ,unsigned INT const maxInt , bool const fStopIt
                 ,const short outCod        , const char *const fName
                 ,void (*rhs)());
  
  int stepperRoss(DOUBLE *y               , DOUBLE *param
               , DOUBLE const x1          , DOUBLE const x2
               , DOUBLE const hInit       , DOUBLE const aTol
               , DOUBLE const rTol        , short const nEdo
               ,unsigned INT const maxInt , bool const fStopIt
               , short const  outCod      , const char *const fName
               , void (*rhs)()            
               , void (*jacY)()           , void (*jacX)());

  int StepperSie(DOUBLE *y                , void **pt
               , DOUBLE const x1          , DOUBLE const x2
               , DOUBLE const hInit       , DOUBLE const hMax
               , DOUBLE const aTol        , DOUBLE const rTol    
               , short const nEdo
               , unsigned INT const maxInt, bool const fStopIt
               , short const  outCod      , const char *const fName
               , void(*rhs)()             , void(*jacY)());

  void writeOutput(int const it     , double const x
                , double const h  , double *y
                , short const nEdo, short const outCod
                , FILE *file );

  double edoError(DOUBLE *y        ,DOUBLE *yOut
               ,DOUBLE *yErr     ,DOUBLE const atol
               ,DOUBLE const rtol,short const nEdo);

  void testAlloc(void* pt,char *const name);

  DOUBLE sqr(DOUBLE const x); 

#endif
/*_EDO_H_*/