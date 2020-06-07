import cantera as ct
import matplotlib.pyplot as plt


def heatC(w,hk):

  return -sum(w*hk)


def CallCantera(fileName,line,t0,stri,nStep,dt=1.e-6,plot='T'):
  gas = ct.Solution(fileName)

  gas.TPY = t0, ct.one_atm, stri

#  r   = ct.IdealGasReactor(gas)
  r   = ct.IdealGasConstPressureReactor(gas)
  sim = ct.ReactorNet([r])
  time  = 0
  fC3H8 = False
  fCH4  = False
  fCO   = False
  fNO   = False
  for name in gas.species_names:
    if name == 'CO':
      fCO = True
      iCO = r.thermo.species_index('CO')
      continue
    if name == 'NO':
      fNO = True
    if name == 'CH4':
      iCH4 = r.thermo.species_index('CH4')
      fCH4 = True
    if name == 'C3H8':
      fC3H8 = True
      iC3H8 = r.thermo.species_index('C3H8')
  
  iO2  = r.thermo.species_index('O2')
  iCO2 = r.thermo.species_index('CO2')
  iH2O = r.thermo.species_index('H2O')


  aO,aH,aC,aN  = r.thermo.atomic_weights
  print('Peso atomico: ')
  print('O : {0}\nH : {1}\n C : {2}\n N: {3}\n'.format(aO,aH,aC,aN))

  mWCH4 = aC + 4*aH
  mWO2  = 2*aO
  mWCO2 =   aC + 2*aO
  mWH2O = 2*aH + aO

  print('Massa Molar: ')
  print('CH4 : {0}\nO2  : {1}\nCO2 : {2}\nH2O: {3}\n'.format(mWCH4,mWO2,mWCO2,mWH2O))

#Arrays to hold the datas
  times = []
  T     = []
  V     = []
  U     = []
  yCH4  = []
  yC3H8 = []
  yO2   = []
  yCO2  = []
  yH2O  = []
  yN2   = []
  yCO   = []
  yNO   = []
  Q     = [0.0]
  q     = 0.0
  QT    = 0.0
  w1    = [0.0]
  w2    = [0.0]
  w3    = [0.0]
  w4    = [0.0]
  H     =[(r.thermo.partial_molar_enthalpies[iCH4]/(1000*mWCH4)\
         ,r.thermo.partial_molar_enthalpies[iO2]/(1000*mWO2)\
         ,r.thermo.partial_molar_enthalpies[iCO2]/(1000*mWCO2)\
         ,r.thermo.partial_molar_enthalpies[iH2O]/(1000*mWH2O))]
  MT    = 0.0
  times.append(time)
# T.append(r.T-273.15)
  T.append(r.T)
  V.append(r.thermo.v)
  U.append(r.thermo.u)
  if fCH4:
    yCH4.append(float(r.thermo['CH4'].Y))
  if fC3H8:
    yC3H8.append(float(r.thermo['C3H8'].Y))
  yO2.append(float(r.thermo['O2'].Y))
  yCO2.append(float(r.thermo['CO2'].Y))
  yH2O.append(float(r.thermo['H2O'].Y))
  yN2.append(float(r.thermo['N2'].Y))
  if fNO:
    yNO. append(float(r.thermo['NO'].Y))
  if fCO:
    yCO.append(float(r.thermo['CO'].Y))

  print('Valores iniciais:')
  print('Pressao   :',r.thermo.P)
  print('cp0(j/kgk):',r.thermo.cp_mass)
  print('cv0(j/kgk):',r.thermo.cv_mass)
  print('Wg        :',r.thermo.mean_molecular_weight)
  print('gama      :',r.thermo.cp_mass/r.thermo.cv_mass)
  print('ro0(kg/m3):',r.thermo.density_mass)
  print('T0 (K)    :',r.thermo.T)

  print('CH4(KJ/Kmol):',r.thermo.partial_molar_enthalpies[iCH4]/1000)
  print('O2 (KJ/Kmol)',r.thermo.partial_molar_enthalpies[iO2]/1000)
  if fCO:
    print('CO (KJ/Kmol)',r.thermo.partial_molar_enthalpies[iCO]/1000)
  print('H2O(KJ/Kmol)',r.thermo.partial_molar_enthalpies[iH2O]/1000)
  print('CO2(KJ/Kmol)',r.thermo.partial_molar_enthalpies[iCO2]/1000)
  
  sim.verbose = True
  for n in range(nStep):
    time+=dt 
    sim.advance(time) 
    times.append(time*1.0e+6)
# temperatura
#    T.append(r.T-273.15)
    T.append(r.T)
# fracao massica
    if fCH4:
      yCH4.append(float(r.thermo['CH4'].Y))
    if fC3H8:
      yC3H8.append(float(r.thermo['C3H8'].Y))
    yO2.append(float(r.thermo['O2'].Y))
    yCO2.append(float(r.thermo['CO2'].Y))
    yH2O.append(float(r.thermo['H2O'].Y))
    yN2.append(float(r.thermo['N2'].Y))
    if fCO:
      yCO.append(float(r.thermo['CO'].Y))
    if fNO:
     yNO.append(float(r.thermo['NO'].Y))
# energia liberada
    q = heatC(r.thermo.net_production_rates
             ,r.thermo.partial_molar_enthalpies)
    Q.append(q/1000) #KJ/s
    QT += q*dt/1000
    MT += r.thermo.net_production_rates[iCH4]*mWCH4*dt
# consumo de ch4 KG/m^3s
    if fCH4:
      w1.append(r.thermo.net_production_rates[iCH4]*mWCH4)
    w2.append(r.thermo.net_production_rates[iO2]*mWO2)
    w3.append(r.thermo.net_production_rates[iCO2]*mWCO2)
    w4.append(r.thermo.net_production_rates[iH2O]*mWH2O)
# entalpia por especie KJ/KG
#    h =(r.thermo.partial_molar_enthalpies[iCH4]/(1000*mWCH4)\
#       ,r.thermo.partial_molar_enthalpies[iO2]/(1000*mWO2)\
#       ,r.thermo.partial_molar_enthalpies[iCO2]/(1000*mWCO2)\
#       ,r.thermo.partial_molar_enthalpies[iH2O]/(1000*mWH2O))
    h =(r.thermo.partial_molar_enthalpies[iCH4]\
       ,r.thermo.partial_molar_enthalpies[iO2]\
       ,r.thermo.partial_molar_enthalpies[iCO2]\
       ,r.thermo.partial_molar_enthalpies[iH2O])
    H.append(h)
    
  print('Valores finais:')
  print('Time      : {0:.8f}'.format(time))  
  print('PressaoF  :',r.thermo.P)
  print('cpF(j/kgk):',r.thermo.cp_mass)
  print('cv0(j/kgk):',r.thermo.cv_mass)
  print('gama      :',r.thermo.cp_mass/r.thermo.cv_mass)
  print('roF(kg/m3):',r.thermo.density_mass)
  print('TF  (K)   :',r.thermo.T)
  print('Q         :',Q[-1])
  print('W1        :',w1[-1])
  print('H0        :',H[0])
  print('H         :',H[-1])
  print('QT(KW/m3) :',QT)
  print('MT        :',-MT)
  if fCH4:
    print('yCH4      :',float(r.thermo['CH4'].Y))
  if fC3H8:
    print('yC3H8      :',float(r.thermo['C3H8'].Y))

  print('Especies com conetracoes maiores que 0.001:')
  for n in gas.species_names:
    if float(r.thermo[n].Y) > 0.0001:  
      print ("{0:10} {1:.4f}".format(n,float(r.thermo[n].Y)))

  plt.xlabel('t(10-6s)')
  if(plot=='T'):
    plt.ylabel('C')
    plt.plot(times, T,lw=1,linestyle=line, label=fileName)
    return
  if(plot=='Q'):
    plt.ylabel('KJ/m^3s')
    plt.plot(times, Q   ,lw=1,linestyle=line,color='black', label=fileName)
    return
  if(plot=='W'):
    plt.ylabel('kg/m^3s')
    plt.semilogx(times,w1,lw=1,linestyle=line,color='black',label='CH4_'+fileName)
#    plt.plot(times,w2,lw=1,linestyle=line,color='blue' ,label='O2_' +fileName)
#    plt.plot(times,w3,lw=1,linestyle=line,color='red'  ,label='CO2_'+fileName)
#    plt.plot(times,w4,lw=1,linestyle=line,color='gray' ,label='CO2_'+fileName)
    return  
  if(plot=='Y'):
    if fCH4:
      plt.plot(times, yCH4,lw=1,linestyle=line,color='sienna', label='yCH4_'+fileName)
    if fC3H8:
      plt.plot(times, yC3H8,lw=1,linestyle=line,color='blue', label='yC3H8_'+fileName)
    plt.plot(times, yO2 ,lw=1,linestyle=line,color='green', label='yO2_'+fileName)
    plt.plot(times, yCO2,lw=1,linestyle=line,color='black', label='yCO2_'+fileName)
    plt.plot(times, yH2O,lw=1,linestyle=line,color='red', label='yH2O_'+fileName)
    plt.plot(times, yN2 ,lw=1,linestyle=line,color='purple', label='yN2_'+fileName)
    if fCO:
      plt.plot(times, yCO ,lw=1,linestyle=line,color='navy', label='yCO_'+fileName)
    if fNO:
      plt.plot(times, yNO ,lw=1,linestyle=line,color='gray', label='yNO_'+fileName)

# ...
    writeRes(times,yCH4,yO2,yH2O,yCO2,yCO,yN2,T)
# .....................................................................................

    return

def writeRes(times,yCH4,yO2,yH2O,yCO2,yCO,yN2,T):

  with open('cantera.out','w') as f:
    for t,CH4,O2,H2O,CO2,CO,N2,T in zip(times,yCH4,yO2,yH2O,yCO2,yCO,yN2,T):
      f.write('{0:.6f} {1:.6f} {2:.6f} {3:.6f} {4:.6f} {5:.6f} {6:.6f} {7:.6f}\n'.format(t,CH4,O2,H2O,CO2,CO,N2,T))


def eq(fileName,stri,t0=298.15):
  gas = ct.Solution(fileName)
  gas.TPY = t0, ct.one_atm, stri
  gas.equilibrate('TP')
  gas()


def steady_state(fileName,T0,stri):
  gas = ct.Solution(fileName)

  gas.TPY = T0, ct.one_atm, stri

  r   = ct.IdealGasReactor(gas)
  sim = ct.ReactorNet([r])
  time  = 0
  fC3H8 = False
  fCH4  = False
  fCO   = False
  fNO   = False
  iCH4  = None
  iC3H8 = None
  for name in gas.species_names:
    if name == 'CO':
      fCO = True
      continue
    if name == 'NO':
      fNO = True
    if name == 'CH4':
      iCH4 = r.thermo.species_index('CH4')
      fCH4 = True
    if name == 'C3H8':
      fC3H8 = True
      iC3H8 = r.thermo.species_index('C3H8')
  
  iO2  = r.thermo.species_index('O2')
  iCO2 = r.thermo.species_index('CO2')
  iH2O = r.thermo.species_index('H2O')
  print(iCH4,iC3H8,iO2,iCO2,iH2O)
#Arrays to hold the datas
  print('Valores iniciais:')
  print('Pressao   :',r.thermo.P)
  print('cp0(j/kgk):',r.thermo.cp_mass)
  print('cv0(j/kgk):',r.thermo.cv_mass)
  print('Wg0       :',r.mean_molecular_weight)
  print('ro0(kg/m3):',r.thermo.density_mass)
  print('T  (K)    :',r.thermo.T-273.15)
  print('H         :',r.thermo.partial_molar_enthalpies*1.e-3)
  sim.verbose = False
  sim.advance_to_steady_state() 
# energia liberada
  q = heatC(r.thermo.net_production_rates
             ,r.thermo.partial_molar_enthalpies)
  Q=q/1000 #KJ/s  
# consumo de ch4
  print('Peso atomico: ')
  print('O : {0}\nH : {1}\n C : {2}\n N: {3}\n'.format(r.thermo.atomic_weights[0]\
                                                    ,r.thermo.atomic_weights[1]\
                                                    ,r.thermo.atomic_weights[2]\
                                                    ,r.thermo.atomic_weights[3]))
    
  print('Valores finais:')
  print('PressaoF  :',r.thermo.P)
  print('cpF(j/kgk):',r.thermo.cp_mass)
  print('cv0(j/kgk):',r.thermo.cv_mass)
  print('Wg0       :',r.mean_molecular_weight)
  print('roF(kg/m3):',r.thermo.density_mass)
  print('QT(KW/m3) :',Q)
  print('T(K)      :',r.thermo.T-273.15)
  print('H         :',r.thermo.partial_molar_enthalpies*1.e-3)

  if fCH4:
    print('yCH4      :',float(r.thermo['CH4'].Y))
  if fC3H8:
    print('yC3H8      :',float(r.thermo['C3H8'].Y))

  print('Especies com conetracoes maiores que 0.001:')
  for n in gas.species_names:
    if float(r.thermo[n].Y) > 0.0001:  
      print ("{0:10} {1:.4f}".format(n,float(r.thermo[n].Y)))


  
CallCantera('ch4_2step_tese.cti','-' ,800.0,'CH4:0.1 , O2:0.6, N2:0.3',100000,dt=1.e-7,plot='Y')
#CallCantera('ch4_cm1.cti','-' ,800.0,'CH4:0.1 , O2:0.6, N2:0.3',100000,dt=1.e-7,plot='Y')

plt.grid()
#plt.ylim([0,-1.e8])
#plt.xlim([1300.0,1500.0])
plt.legend()
plt.show()

#eq('c3h8_cm1.cti','C3H8:0.005, O2:0.295, N2:0.7')

#steady_state('sandiego.cti'   ,295.15,'C3H8:0.1, O2:0.2, N2:0.7')
#steady_state('ch4_cm1.cti'   ,298.15,'CH4:0.055, O2:0.22, N2:0.725')
#steady_state('ch4_wd.cti'    ,298.15,'CH4:0.055, O2:0.22, N2:0.725')



