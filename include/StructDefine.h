#ifndef _STRUCTDEFINE_H_
  #define _STRUCTDEFINE_H_
  
  #include<Define.h>
  #include<NasaPol.h>


  typedef struct 
  {
    DOUBLE A,E,Ta,beta;  
  
  } Arrhenius;

  typedef struct 
  {
    PolNasa pol[MAX_SPECIES];  
  
  } PropVar;

  typedef struct {
  
    bool reverse;
  
    DOUBLE stch[MAX_SPECIES][2];
    DOUBLE exp[MAX_SPECIES][2];
  
    Arrhenius ArrF,ArrR; 

  } Reaction;

  typedef struct 
  {
    char  name[MAX_SPECIES][MAX_NAME_SP];
    unsigned nReac,nSp; 
    DOUBLE mW[MAX_SPECIES];
    Reaction reac[MAX_REAC]; 
  } Chemical;

#endif
