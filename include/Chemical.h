#ifndef _CHEMICAL_H_
  #define _CHEMICAL_H_

  #include<math.h>
  #include<Define.h>  
  #include<StructDefine.h>

  void chemicalRsh(DOUBLE const t    ,DOUBLE *RESTRICT y
                ,DOUBLE *RESTRICT w,void **pt);

  DOUBLE yFracSum(DOUBLE *RESTRICT yFrac, unsigned short n) ;

  DOUBLE molarMass(DOUBLE *RESTRICT y     ,DOUBLE *RESTRICT mW
                  ,unsigned short const nSp);
#endif /*_CHEMICAL_H_*/
