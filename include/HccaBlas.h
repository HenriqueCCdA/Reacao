#ifndef _HCCABLAS_
  #define _HCCABLAS_
  #define NFUNC           159
  #define HCCABLASZERO 1.0e-8
  #include<Define.h>
  #include<stdio.h>
  #ifdef _OPENMP
    #include<omp.h>
    DOUBLE tmpDotOmp;
    DOUBLE tmpDotOmp1;
    DOUBLE tmpDotOmp2;
    DOUBLE tmpDotOmp3;
    DOUBLE tmpDotOmp4;
    DOUBLE tmpDotOmp5;
    DOUBLE tmpDotOmp6;
    DOUBLE tmpDotOmp7;
    DOUBLE tmpDotOmp8;
  #endif

/*...*/
  void getNameHccaBlas(void);
  long  lopMatVecFull(INT nLin,INT nCol);
  long flopMatVecCsr(INT neq,INT nad,short ty);
  long flopDot(INT nDim);
/*produto vetorial*/  
  void prodVet(DOUBLE *RESTRICT a,DOUBLE *RESTRICT b
              ,DOUBLE *RESTRICT c);
  INT xDiffY(DOUBLE *RESTRICT x,DOUBLE *RESTRICT y
            ,DOUBLE const tol  , INT n);
/*...................................................................*/

/*... normas*/
  LDOUBLE lnormInf(LDOUBLE *RESTRICT x,INT const lin,INT const col);

/*...................................................................*/

/*level 1*/
  void alphaProdVector(DOUBLE const alpha,DOUBLE *RESTRICT a
                      ,INT const nDim    ,DOUBLE *RESTRICT c); 
  void addVector(DOUBLE const alpha,DOUBLE *RESTRICT a
                ,DOUBLE const beta ,DOUBLE *RESTRICT b
                ,INT const nDim    ,DOUBLE *RESTRICT c); 

  DOUBLE dot(DOUBLE *RESTRICT x1,DOUBLE *RESTRICT x2,INT n);
  DOUBLE dotO2L2(DOUBLE *RESTRICT x1,DOUBLE *RESTRICT x2,INT n);
  DOUBLE   dotL2(DOUBLE *RESTRICT x1,DOUBLE *RESTRICT x2,INT n);
  DOUBLE   dotL4(DOUBLE *RESTRICT x1,DOUBLE *RESTRICT x2,INT n);
  DOUBLE   dotL6(DOUBLE *RESTRICT x1,DOUBLE *RESTRICT x2,INT n);
  DOUBLE   dotL8(DOUBLE *RESTRICT x1,DOUBLE *RESTRICT x2,INT n);
  DOUBLE   dotO2(DOUBLE *RESTRICT x1,DOUBLE *RESTRICT x2,INT n);
  DOUBLE   dotO4(DOUBLE *RESTRICT x1,DOUBLE *RESTRICT x2,INT n);
  DOUBLE   dotO6(DOUBLE *RESTRICT x1,DOUBLE *RESTRICT x2,INT n);
  DOUBLE   dotO8(DOUBLE *RESTRICT x1,DOUBLE *RESTRICT x2,INT n);
#if _OPENMP
  DOUBLE     dotOmp(DOUBLE *RESTRICT x1,DOUBLE *RESTRICT x2,INT n);
  DOUBLE dotOmpO2L2(DOUBLE *RESTRICT x1,DOUBLE *RESTRICT x2,INT n);
  DOUBLE dotOmpO2L4(DOUBLE *RESTRICT x1,DOUBLE *RESTRICT x2,INT n);
  DOUBLE   dotOmpL2(DOUBLE *RESTRICT x1,DOUBLE *RESTRICT x2,INT n);
  DOUBLE   dotOmpL4(DOUBLE *RESTRICT x1,DOUBLE *RESTRICT x2,INT n);
  DOUBLE   dotOmpL6(DOUBLE *RESTRICT x1,DOUBLE *RESTRICT x2,INT n);
  DOUBLE   dotOmpL8(DOUBLE *RESTRICT x1,DOUBLE *RESTRICT x2,INT n);
  DOUBLE   dotOmpO2(DOUBLE *RESTRICT x1,DOUBLE *RESTRICT x2,INT n);
  DOUBLE   dotOmpO4(DOUBLE *RESTRICT x1,DOUBLE *RESTRICT x2,INT n);
  DOUBLE   dotOmpO6(DOUBLE *RESTRICT x1,DOUBLE *RESTRICT x2,INT n);
  DOUBLE   dotOmpO8(DOUBLE *RESTRICT x1,DOUBLE *RESTRICT x2,INT n);
#endif
/*...................................................................*/

/*level 2*/
/* ... matriz cheia*/
  void matVecFull(DOUBLE *RESTRICT a
                 ,DOUBLE *RESTRICT x
                 ,DOUBLE *RESTRICT y
                 ,INT nLin          ,INT nCol);
  void matVecFullO2(DOUBLE *RESTRICT a
                 ,DOUBLE *RESTRICT x
                 ,DOUBLE *RESTRICT y
                 ,INT nLin          ,INT nCol);
  void matVecFullO4(DOUBLE *RESTRICT a
                   ,DOUBLE *RESTRICT x
                   ,DOUBLE *RESTRICT y
                   ,INT nLin          ,INT nCol);
  void matVecFullO2I2(DOUBLE *RESTRICT a
                     ,DOUBLE *RESTRICT x
                     ,DOUBLE *RESTRICT y
                     ,INT nLin          ,INT nCol);
  void matVecFullO4I2(DOUBLE *RESTRICT a
                     ,DOUBLE *RESTRICT x
                     ,DOUBLE *RESTRICT y
                     ,INT nLin          ,INT nCol);
  void matVecFullO4I4(DOUBLE *RESTRICT a
                     ,DOUBLE *RESTRICT x
                     ,DOUBLE *RESTRICT y
                     ,INT nLin          ,INT nCol);
  void lmatVecFull(LDOUBLE *RESTRICT a
                  ,LDOUBLE *RESTRICT x
                  ,LDOUBLE *RESTRICT y
                  ,INT nLin       
                  ,INT nCol);
/*...................................................................*/

/*...Csr*/
//typedef enum {csr=1,csrD=2,csrC=3} typeCsr;

/*... CsrD*/ 
  void     matVecCsrD(INT neq           
                     ,INT *RESTRICT ia  ,INT *RESTRICT ja
                     ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT ad
                     ,DOUBLE *RESTRICT x,DOUBLE *RESTRICT y);
  void   matVecCsrDI2(INT neq           
                     ,INT *RESTRICT ia  ,INT *RESTRICT ja
                     ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT ad
                     ,DOUBLE *RESTRICT x,DOUBLE *RESTRICT y);
  void   matVecCsrDI4(INT neq           
                     ,INT *RESTRICT ia  ,INT *RESTRICT ja
                     ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT ad
                     ,DOUBLE *RESTRICT x,DOUBLE *RESTRICT y);
  void   matVecCsrDI6(INT neq           
                     ,INT *RESTRICT ia  ,INT *RESTRICT ja
                     ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT ad
                     ,DOUBLE *RESTRICT x,DOUBLE *RESTRICT y);
  void   matVecCsrDO2(INT neq           
                     ,INT *RESTRICT ia  ,INT *RESTRICT ja
                     ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT ad
                     ,DOUBLE *RESTRICT x,DOUBLE *RESTRICT y);
  void   matVecCsrDO4(INT neq           
                     ,INT *RESTRICT ia  ,INT *RESTRICT ja
                     ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT ad
                     ,DOUBLE *RESTRICT x,DOUBLE *RESTRICT y);
  void   matVecCsrDO6(INT neq           
                     ,INT *RESTRICT ia  ,INT *RESTRICT ja
                     ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT ad
                     ,DOUBLE *RESTRICT x,DOUBLE *RESTRICT y);
  void matVecCsrDO2I2(INT neq           
                     ,INT *RESTRICT ia  ,INT *RESTRICT ja
                     ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT ad
                     ,DOUBLE *RESTRICT x,DOUBLE *RESTRICT y);
#ifdef _OPENMP
/*... CsrDOmp*/ 
  void       matVecCsrDOmp(INT neq           
                     ,INT *RESTRICT ia  ,INT *RESTRICT ja
                     ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT ad
                     ,DOUBLE *RESTRICT x,DOUBLE *RESTRICT y);
  void       matVecCsrDOmpI2(INT neq           
                     ,INT *RESTRICT ia  ,INT *RESTRICT ja
                     ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT ad
                     ,DOUBLE *RESTRICT x,DOUBLE *RESTRICT y);
  void       matVecCsrDOmpI4(INT neq           
                     ,INT *RESTRICT ia  ,INT *RESTRICT ja
                     ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT ad
                     ,DOUBLE *RESTRICT x,DOUBLE *RESTRICT y);
  void       matVecCsrDOmpI6(INT neq           
                     ,INT *RESTRICT ia  ,INT *RESTRICT ja
                     ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT ad
                     ,DOUBLE *RESTRICT x,DOUBLE *RESTRICT y);
  void       matVecCsrDOmpO2(INT neq           
                     ,INT *RESTRICT ia  ,INT *RESTRICT ja
                     ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT ad
                     ,DOUBLE *RESTRICT x,DOUBLE *RESTRICT y);
  void       matVecCsrDOmpO4(INT neq           
                     ,INT *RESTRICT ia  ,INT *RESTRICT ja
                     ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT ad
                     ,DOUBLE *RESTRICT x,DOUBLE *RESTRICT y);
  void       matVecCsrDOmpO6(INT neq           
                     ,INT *RESTRICT ia  ,INT *RESTRICT ja
                     ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT ad
                     ,DOUBLE *RESTRICT x,DOUBLE *RESTRICT y);
  void       matVecCsrDOmpO2I2(INT neq           
                     ,INT *RESTRICT ia  ,INT *RESTRICT ja
                     ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT ad
                     ,DOUBLE *RESTRICT x,DOUBLE *RESTRICT y);
/*balanciamento manual do trabalho entre as threads*/  
  void matVecCsrDOmpBal(INT neq           
                     ,INT *RESTRICT ia  ,INT *RESTRICT ja
                     ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT ad
                     ,DOUBLE *RESTRICT x,DOUBLE *RESTRICT y
                     ,INT  *thBegin     ,INT *thEnd  );
  void matVecCsrDOmpBalI2(INT neq           
                     ,INT *RESTRICT ia  ,INT *RESTRICT ja
                     ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT ad
                     ,DOUBLE *RESTRICT x,DOUBLE *RESTRICT y
                     ,INT  *thBegin     ,INT *thEnd  );
  void matVecCsrDOmpBalI4(INT neq           
                     ,INT *RESTRICT ia  ,INT *RESTRICT ja
                     ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT ad
                     ,DOUBLE *RESTRICT x,DOUBLE *RESTRICT y
                     ,INT  *thBegin     ,INT *thEnd  );
  void matVecCsrDOmpBalI6(INT neq           
                     ,INT *RESTRICT ia  ,INT *RESTRICT ja
                     ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT ad
                     ,DOUBLE *RESTRICT x,DOUBLE *RESTRICT y
                     ,INT  *thBegin     ,INT *thEnd  );
  void matVecCsrDOmpBalO2(INT neq           
                     ,INT *RESTRICT ia  ,INT *RESTRICT ja
                     ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT ad
                     ,DOUBLE *RESTRICT x,DOUBLE *RESTRICT y
                     ,INT  *thBegin     ,INT *thEnd  );
  void matVecCsrDOmpBalO4(INT neq           
                     ,INT *RESTRICT ia  ,INT *RESTRICT ja
                     ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT ad
                     ,DOUBLE *RESTRICT x,DOUBLE *RESTRICT y
                     ,INT  *thBegin     ,INT *thEnd  );
  void matVecCsrDOmpBalO6(INT neq           
                     ,INT *RESTRICT ia  ,INT *RESTRICT ja
                     ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT ad
                     ,DOUBLE *RESTRICT x,DOUBLE *RESTRICT y
                     ,INT  *thBegin     ,INT *thEnd  );
  void matVecCsrDOmpBalO2I2(INT neq           
                     ,INT *RESTRICT ia  ,INT *RESTRICT ja
                     ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT ad
                     ,DOUBLE *RESTRICT x,DOUBLE *RESTRICT y
                     ,INT  *thBegin     ,INT *thEnd  );


#endif         
/*... CsrC*/ 
  void     matVecCsrC(INT neq           
                     ,INT *RESTRICT ia   ,INT *RESTRICT ja
                     ,DOUBLE *RESTRICT au,DOUBLE *RESTRICT ad
                     ,DOUBLE *RESTRICT al
                     ,DOUBLE *RESTRICT x ,DOUBLE *RESTRICT y);
  void   matVecCsrCI2(INT neq           
                     ,INT *RESTRICT ia   ,INT *RESTRICT ja
                     ,DOUBLE *RESTRICT au,DOUBLE *RESTRICT ad
                     ,DOUBLE *RESTRICT al
                     ,DOUBLE *RESTRICT x ,DOUBLE *RESTRICT y);
  void   matVecCsrCI4(INT neq           
                     ,INT *RESTRICT ia   ,INT *RESTRICT ja
                     ,DOUBLE *RESTRICT au,DOUBLE *RESTRICT ad
                     ,DOUBLE *RESTRICT al
                     ,DOUBLE *RESTRICT x ,DOUBLE *RESTRICT y);
  void   matVecCsrCI6(INT neq           
                     ,INT *RESTRICT ia   ,INT *RESTRICT ja
                     ,DOUBLE *RESTRICT au,DOUBLE *RESTRICT ad
                     ,DOUBLE *RESTRICT al
                     ,DOUBLE *RESTRICT x ,DOUBLE *RESTRICT y);
  void   matVecCsrCO2(INT neq           
                     ,INT *RESTRICT ia   ,INT *RESTRICT ja
                     ,DOUBLE *RESTRICT au,DOUBLE *RESTRICT ad
                     ,DOUBLE *RESTRICT al
                     ,DOUBLE *RESTRICT x ,DOUBLE *RESTRICT y);
  void   matVecCsrCO4(INT neq           
                     ,INT *RESTRICT ia   ,INT *RESTRICT ja
                     ,DOUBLE *RESTRICT au,DOUBLE *RESTRICT ad
                     ,DOUBLE *RESTRICT al
                     ,DOUBLE *RESTRICT x ,DOUBLE *RESTRICT y);
  void   matVecCsrCO6(INT neq           
                     ,INT *RESTRICT ia   ,INT *RESTRICT ja
                     ,DOUBLE *RESTRICT au,DOUBLE *RESTRICT ad
                     ,DOUBLE *RESTRICT al
                     ,DOUBLE *RESTRICT x ,DOUBLE *RESTRICT y);
  void matVecCsrCO2I2(INT neq           
                     ,INT *RESTRICT ia   ,INT *RESTRICT ja
                     ,DOUBLE *RESTRICT au,DOUBLE *RESTRICT ad
                     ,DOUBLE *RESTRICT al
                     ,DOUBLE *RESTRICT x ,DOUBLE *RESTRICT y);
#ifdef _OPENMP
  void     matVecCsrCOmp(INT neq           
                        ,INT *RESTRICT ia    ,INT *RESTRICT ja
                        ,DOUBLE *RESTRICT au ,DOUBLE *RESTRICT ad
                        ,DOUBLE *RESTRICT al
                        ,DOUBLE *RESTRICT x  ,DOUBLE *RESTRICT y
                        ,INT  *thBegin       ,INT *thEnd  
                        ,INT  *thHeight    
                        ,DOUBLE *RESTRICT thY,int nThreads);          
  void     matVecCsrCOmpI2(INT neq           
                          ,INT *RESTRICT ia    ,INT *RESTRICT ja
                          ,DOUBLE *RESTRICT au ,DOUBLE *RESTRICT ad
                          ,DOUBLE *RESTRICT al
                          ,DOUBLE *RESTRICT x  ,DOUBLE *RESTRICT y
                          ,INT  *thBegin       ,INT *thEnd  
                          ,INT  *thHeight    
                          ,DOUBLE *RESTRICT thY,int nThreads);          
  void     matVecCsrCOmpI4(INT neq           
                          ,INT *RESTRICT ia    ,INT *RESTRICT ja
                          ,DOUBLE *RESTRICT au ,DOUBLE *RESTRICT ad
                          ,DOUBLE *RESTRICT al
                          ,DOUBLE *RESTRICT x  ,DOUBLE *RESTRICT y
                          ,INT  *thBegin       ,INT *thEnd  
                          ,INT  *thHeight    
                          ,DOUBLE *RESTRICT thY,int nThreads);          
  void     matVecCsrCOmpI6(INT neq           
                          ,INT *RESTRICT ia    ,INT *RESTRICT ja
                          ,DOUBLE *RESTRICT au ,DOUBLE *RESTRICT ad
                          ,DOUBLE *RESTRICT al
                          ,DOUBLE *RESTRICT x  ,DOUBLE *RESTRICT y
                          ,INT  *thBegin       ,INT *thEnd  
                          ,INT  *thHeight    
                          ,DOUBLE *RESTRICT thY,int nThreads);          
  void     matVecCsrCOmpO2(INT neq           
                          ,INT *RESTRICT ia    ,INT *RESTRICT ja
                          ,DOUBLE *RESTRICT au ,DOUBLE *RESTRICT ad
                          ,DOUBLE *RESTRICT al
                          ,DOUBLE *RESTRICT x  ,DOUBLE *RESTRICT y
                          ,INT  *thBegin       ,INT *thEnd  
                          ,INT  *thHeight    
                          ,DOUBLE *RESTRICT thY,int nThreads);          
  void     matVecCsrCOmpO4(INT neq           
                          ,INT *RESTRICT ia    ,INT *RESTRICT ja
                          ,DOUBLE *RESTRICT au ,DOUBLE *RESTRICT ad
                          ,DOUBLE *RESTRICT al
                          ,DOUBLE *RESTRICT x  ,DOUBLE *RESTRICT y
                          ,INT  *thBegin       ,INT *thEnd  
                          ,INT  *thHeight    
                          ,DOUBLE *RESTRICT thY,int nThreads);          
  void     matVecCsrCOmpO6(INT neq           
                          ,INT *RESTRICT ia    ,INT *RESTRICT ja
                          ,DOUBLE *RESTRICT au ,DOUBLE *RESTRICT ad
                          ,DOUBLE *RESTRICT al
                          ,DOUBLE *RESTRICT x  ,DOUBLE *RESTRICT y
                          ,INT  *thBegin       ,INT *thEnd  
                          ,INT  *thHeight    
                          ,DOUBLE *RESTRICT thY,int nThreads);          
  void     matVecCsrCOmpO2I2(INT neq           
                            ,INT *RESTRICT ia    ,INT *RESTRICT ja
                            ,DOUBLE *RESTRICT au ,DOUBLE *RESTRICT ad
                            ,DOUBLE *RESTRICT al
                            ,DOUBLE *RESTRICT x  ,DOUBLE *RESTRICT y
                            ,INT  *thBegin       ,INT *thEnd  
                            ,INT  *thHeight    
                            ,DOUBLE *RESTRICT thY,int nThreads);          
#endif         
/*... Csr*/ 
  void     matVecCsr(INT neq           
                    ,INT *RESTRICT ia  ,INT *RESTRICT ja
                    ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT x
                    ,DOUBLE *RESTRICT y);
  void   matVecCsrI2(INT neq           
                    ,INT *RESTRICT ia  ,INT *RESTRICT ja
                    ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT x
                    ,DOUBLE *RESTRICT y);
  void   matVecCsrI4(INT neq           
                    ,INT *RESTRICT ia  ,INT *RESTRICT ja
                    ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT x
                    ,DOUBLE *RESTRICT y);
  void   matVecCsrI6(INT neq           
                    ,INT *RESTRICT ia  ,INT *RESTRICT ja
                    ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT x
                    ,DOUBLE *RESTRICT y);
  void   matVecCsrO2(INT neq           
                    ,INT *RESTRICT ia  ,INT *RESTRICT ja
                    ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT x
                    ,DOUBLE *RESTRICT y);
  void   matVecCsrO4(INT neq           
                    ,INT *RESTRICT ia  ,INT *RESTRICT ja
                    ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT x
                    ,DOUBLE *RESTRICT y);
  void   matVecCsrO6(INT neq           
                    ,INT *RESTRICT ia  ,INT *RESTRICT ja
                    ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT x
                    ,DOUBLE *RESTRICT y);
  void matVecCsrO2I2(INT neq           
                    ,INT *RESTRICT ia  ,INT *RESTRICT ja
                    ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT x
                    ,DOUBLE *RESTRICT y);
#ifdef _OPENMP
/*... CsrOmp*/ 
  void       matVecCsrOmp(INT neq           
                         ,INT *RESTRICT ia  ,INT *RESTRICT ja
                         ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT x
                         ,DOUBLE *RESTRICT y);
  void       matVecCsrOmpI2(INT neq           
                           ,INT *RESTRICT ia  ,INT *RESTRICT ja
                           ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT x
                           ,DOUBLE *RESTRICT y);
  void       matVecCsrOmpI4(INT neq           
                           ,INT *RESTRICT ia  ,INT *RESTRICT ja
                           ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT x
                           ,DOUBLE *RESTRICT y);
  void       matVecCsrOmpI6(INT neq           
                           ,INT *RESTRICT ia  ,INT *RESTRICT ja
                           ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT x
                           ,DOUBLE *RESTRICT y);
  void       matVecCsrOmpO2(INT neq           
                           ,INT *RESTRICT ia  ,INT *RESTRICT ja
                           ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT x
                           ,DOUBLE *RESTRICT y);
  void       matVecCsrOmpO4(INT neq           
                           ,INT *RESTRICT ia  ,INT *RESTRICT ja
                           ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT x
                           ,DOUBLE *RESTRICT y);
  void       matVecCsrOmpO6(INT neq           
                           ,INT *RESTRICT ia  ,INT *RESTRICT ja
                           ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT x
                           ,DOUBLE *RESTRICT y);
  void     matVecCsrOmpO2I2(INT neq           
                           ,INT *RESTRICT ia  ,INT *RESTRICT ja
                           ,DOUBLE *RESTRICT a,DOUBLE *RESTRICT x
                           ,DOUBLE *RESTRICT y);
#endif
/*...................................................................*/
#endif
