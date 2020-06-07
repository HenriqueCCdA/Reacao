import matplotlib.pyplot as plt

with open('teste.out','r') as f:
    data = f.read()
        
lines = data.splitlines()[:]

step = []
x  = []
h  = []
y1 = []
y2 = []
y3 = []
y4 = []
y5 = []
y6 = []
T  = [] 
for line in lines:
    word = line.split()
    step.append(int(word[0]))
    x.append(float(word[1]))
    h.append(float(word[2]))
    y1.append(float(word[3]))
    y2.append(float(word[4]))
    y3.append(float(word[5]))
    y4.append(float(word[6]))
    y5.append(float(word[7]))
    y6.append(float(word[8]))
#    T.append(float(word[9])-273.15)
    T.append(float(word[9]))

with open('cantera.out','r') as f:
    data = f.read()
        
lines = data.splitlines()[:]

cx  = []
cy1 = []
cy2 = []
cy3 = []
cy4 = []
cy5 = []
cy6 = []
cT  = [] 
for line in lines:
    word = line.split()
    cx.append(float(word[0])/1.e+06)
    cy1.append(float(word[1]))
    cy2.append(float(word[2]))
    cy3.append(float(word[3]))
    cy4.append(float(word[4]))
    cy5.append(float(word[5]))
    cy6.append(float(word[6]))
    cT.append(float(word[7]))
#
fig, ax1 = plt.subplots(figsize=(10,10))

ax1.set_xlabel('t(s)')

ax2 = ax1.twinx()
ax2.set_ylabel('Temperature(K)')
ax1.set_ylabel('Y')


ax2.set_ylim([500,3500])
plt.xlim([0.001,0.004])
#plt.xlim([0,100])
ax1.plot(x ,y1 ,linestyle='-',label='CH4-IE'     ,color='sienna')
ax1.plot(cx,cy1,linestyle='--',label='CH4-Cantera',color='sienna', marker='^', markevery=1000)
#
ax1.plot(x ,y2 ,linestyle='-',label='O2-IE     ',color='green')
ax1.plot(cx,cy2,linestyle='--',label='O2-Cantera',color='darkgreen', marker='s', markevery=1000)
#
ax1.plot(x ,y3 ,linestyle='-',label='H2O-IE     ',color='red')
ax1.plot(cx,cy3,linestyle='--',label='H2O-Cantera',color='darkred', marker='v', markevery=1000)
#
ax1.plot(x ,y4 ,linestyle='-',label='CO2-IE     ',color='gray')
ax1.plot(cx,cy4,linestyle='--',label='CO2-Cantera',color='gray', marker='o', markevery=1000)
#
ax1.plot(x ,y5 ,linestyle='-',label='CO-IE     ',color='m')
ax1.plot(cx,cy5,linestyle='--',label='CO-Cantera',color='m', marker='x', markevery=1000)
#
ax1.plot(x ,y6 ,linestyle='-',label='N2-IE     ',color='blue')
ax1.plot(cx,cy6,linestyle=':',label='N2-Cantera',color='darkblue', marker='p', markevery=1000)
    
#
ax2.plot(x  ,T,linestyle='-',label='Temperature-IE     ',color='black')
ax2.plot(cx,cT,linestyle='--',label='Temperature-Cantera',color='black', marker='h', markevery=1000)

ax1.legend(bbox_to_anchor=(0.8, 0.55))
ax2.legend(bbox_to_anchor=(0.3, 0.9))
#plt.show()
plt.savefig('ch4_2step.pdf')
#plt.plot(x,h,label='h')
#plt.show()
