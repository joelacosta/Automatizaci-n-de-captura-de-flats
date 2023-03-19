#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt

#Defino los parametros del modelo. La mayoria son parametros estructurales del puente. Sigma es la desviacion de la distruibucion de frencuencias con las que peatos caminan.
K=4.73e6
M=1.13e5
B=1.1e4
w_0=np.sqrt(K/M)
G=30
C=1.05
sigma=0.63

#Esto define el metodo de runge-kutta de orden 4
def rk4(x, dt):
    x = np.asarray(x)
    k1 = campo_vector(x)*dt
    k2 = campo_vector(x + k1*0.5)*dt
    k3 = campo_vector(x + k2*0.5)*dt
    k4 = campo_vector(x + k3)*dt
    return x + (k1 + 2*k2 + 2*k3 + k4)/6

#Aca defino el campo vector. Tiene n+2 lugares, el primero para X el segundo para X' y los dempas para los n peatones.
def campo_vector(z):
    n=len(z)-2
    d=np.zeros(len(z))
    suma =0
    for i in range(n):
        suma=suma+np.sin(z[i+2])
    d[0] = z[1]
    d[1] = -B*z[1]/M - K*z[0]/M+G/M*suma
    for i in range(n):
        d[i+2]=omega[i]-C*np.sin(z[i+2])*d[1]
    return d

#Dado un vector de condiciones iniciales x_0 en t_0, integra las ecuaciones hasta t_0+delta_t para cantidad de peatones igual a la que habia antes +delta_n (delta_n puede ser 0). Tomas las frecuencias w_i y el paso en el tiempo que esten definidas globalmente
def solucion(t_0,delta_t,x_0,delta_n):
    t= np.arange(t_0, t_0+delta_t, step=dt)
    x_0=x_0 + (2*np.pi*np.random.random(delta_n)).tolist()
    sol=[0 for i in t]
    sol[0]=x_0
    for i in range(len(t)-1):
        sol[i+1]=rk4(sol[i],dt).tolist()
    return t.tolist(),sol

#Le das las soluciones que salen de la funcion anterior, y te devuelve una listas con la cantidad de peatones en funcion del tiempo, la amplitud y la sincronización
def Amp_Fase(sol):
    A=[0 for i in sol]
    fase=[0 for i in sol]
    N=[0 for i in sol]
    for i in range(len(sol)):
        A[i]=np.sqrt(sol[i][0]**2+sol[i][1]**2/w_0**2)
        k=len(sol[i])
        N[i]=k-2
        s=0
        for j in range(N[i]):
            s=s+np.exp(1j*sol[i][j+2])
        fase[i]=np.abs(s)/(k-1)
    return N,A,fase

#Esto randomiza las frecuencias con las que cada peaton camina siguiendo una distribucion normal. El valor de sigma esta muy estudiado.
def omegas(n,sigma):
    return np.random.normal(w_0,sigma,n)

#Ejemplo de solucion. Toma 300 peatones en el punte y evoluciona el sistema desde t_0=0 hasta t_f=300
n=300
omega=omegas(n,sigma)
x_0=[0,0]
t,sol=solucion(0,300,x_0,n)
plt.plot(t,[i[0] for i in sol])


#Partiendo de n=50 resuelvo sumando en cada paso 10 peatones y usando las condiciones finales del paso anterior como inicales de este.
#Esto era para simular la figura 2 del paper de Strogatz
n=300
omega=omegas(n,sigma)
x_0=[0,0]+(np.random.random(50)).tolist()
sol=[x_0]
t=[0]
delta_n=10
dt=0.1
for i in range(17):
    t_i,sol_i=solucion(t[-1],150,sol[-1],delta_n)
    t=t+t_i
    sol=sol+sol_i
    print(i*delta_n) #para ver por donde va
N_2,Amp_2,Fase_2=Amp_Fase(sol)
t_2=t

#Grafico lo anterior
fig, ax = plt.subplots(3,figsize=(10,12),sharex=True)
label_size=16
title_size=18

tc=t_2[N_2.index(150)]
ax[0].plot(t_2,N_2)
ax[0].grid('on')
ax[0].set_ylabel('Cantidad de\n peatones',fontsize=label_size)
ax[0].set_ylim([50,N_2[-1]+10])
ax[0].axvline(tc,color='k',linestyle='dashed')

ax[1].plot(t_2,[i*100 for i in Amp_2])
ax[1].grid('on')
ax[1].set_ylabel('Amplitud del puente\n (cm)', fontsize=label_size)
ax[1].set_ylim([-0.5,100*Amp_2[-1]+0.5])
ax[1].axvline(tc,color='k',linestyle='dashed')

ax[2].plot(t_2,Fase_2)
ax[2].grid('on')
ax[2].set_xlabel('Tiempo [s]',fontsize=label_size)
ax[2].set_ylabel(u'Sincronización entre \n peatones',fontsize=label_size)
ax[2].set_ylim([0,1])
ax[2].axvline(tc,color='k',linestyle='dashed')
plt.xlim([0,t_2[-1]])
plt.savefig('fig2.png',dpi=300,bbox_inches='tight',facecolor='white')


#Simulo para 300
n=300
omega=(np.sort(omegas(n,sigma))).tolist()
dt=0.03
sol=[[0,0]]
t=[0]
t_d,sol_d=solucion(t[-1],100,sol[0],n)
N_d,Amp_d,Fase_d=Amp_Fase(sol_d)

#Grafico
fases=[i[2:] for i in sol_d]
fases_d=(np.cos(np.array(fases))).transpose()
x_d,y_d=np.meshgrid(t_d,[i for i in range(n)])

label_size=18
f, ax = plt.subplots(2, gridspec_kw={'height_ratios': [1, 3]},figsize=(10,10),sharex=True)
ax[0].plot(t_d,[100*i[0] for i in sol_d])
ax[0].set_ylim([-10,10])
ax[0].set_ylabel('Amplitud \n(cm)',fontsize=label_size)
ax[0].grid('on')

ax[1].pcolormesh(x_d,y_d,fases_d[:-1,:-1],cmap='Blues')
ax[1].set_ylabel('Índice de cada peatón',fontsize=label_size)
n_0=[omega.index(i) for i in omega if i>w_0][0]
ax[1].axhline(n_0,color='k',alpha=0.7)
plt.xlabel('Tiempo (s)', fontsize=label_size)
plt.xlim([0,60])
plt.savefig('n300.png',dpi=600,bbox_inches='tight',facecolor='white')

#Simulo para 100
n=100
omega=np.random.normal(w_0,sigma,n)
dt=0.05
sol=[[0,0]]
t=[0]
t_u,sol_u=solucion(t[-1],200,sol[0],n)
N_u,Amp_u,Fase_u=Amp_Fase(sol_u)


#Grafico
fases=[i[2:] for i in sol_u]
fases_u=(np.cos(np.array(fases))).transpose()
x_u,y_u=np.meshgrid(t_u,[i for i in range(n)])
label_size=18

f, ax = plt.subplots(2, gridspec_kw={'height_ratios': [1, 3]},figsize=(10,10),sharex=True)
ax[0].plot(t_u,[100*i[0] for i in sol_u])
ax[0].set_ylim([-10,10])
ax[0].set_ylabel('Amplitud \n(cm)',fontsize=label_size)
ax[0].grid('on')

ax[1].pcolormesh(x_u,y_u,fases_u[:-1,:-1],cmap='Blues')
ax[1].set_ylabel('Índice de cada peatón',fontsize=label_size)
plt.xlabel('Tiempo (s)', fontsize=label_size)
plt.xlim([0,60])
plt.savefig('n130.png',dpi=600,bbox_inches='tight',facecolor='white')


#Comparo los resultados entre los dos anteriores
f,ax=plt.subplots(2,figsize=(8,8),sharex=True)
label_size=14
ax[0].plot(t_u,[100*i for i in Amp_u],label='N = 130')
ax[0].plot(t_d,[100*i for i in Amp_d],label='N = 300')
ax[0].plot(t_d,[100*i[0] for i in sol_d],color='orange',alpha=0.5,zorder=3)
ax[0].plot(t_u,[100*i[0] for i in sol_u],color='blue',zorder=4)
ax[0].legend()
ax[0].grid('on')
ax[0].set_ylabel('Amplitud del puente\n (cm)', fontsize=label_size)
ax[0].set_ylim([-1,12])

ax[1].plot(t_u,Fase_u)
ax[1].plot(t_d,Fase_d)
ax[1].grid('on')
ax[1].set_ylabel(u'Sincronización entre \n peatones',fontsize=label_size)
plt.xlim([0,100])
plt.xlabel('Tiempo (s)',fontsize=label_size)
plt.savefig('comparacion.png',dpi=300,bbox_inches='tight',facecolor='white')


f=plt.figure(figsize=(8,6))
ax=f.add_subplot(111)
for j in range(n):
    if n%30==0:
        fase=[i[j+2] for i in sol_d]
        #plt.plot(t_d,np.asarray(fase) -[w_0*i for i in t_d])
        plt.plot(t_d,(np.asarray(fase) +tita)/np.pi)


plt.grid('on')
plt.yticks([(i) for i in range(-12,13)],[str(i)+r'$\pi$' for i in range(-12,13)])
ax.yaxis.tick_right()
ax.yaxis.set_label_position("right")
plt.xlim([0,80])
plt.ylim([-10,10])
plt.ylabel(r'$\Theta_i - \psi$',fontsize=label_size)
plt.xlabel('Tiempo (s)',fontsize=label_size)
plt.savefig('fases.png',dpi=300,bbox_inches='tight',facecolor='white')

#Aca comparo las soluciones con una solucion clasica del tipo sen. 
fuerza=[np.sin(i[150]) for i in sol_d]
t_f=t_d[-500:]
fuerza_f=fuerza[-500:]

x_f=[i[0] for i in sol_d[-500:]]
x_max=np.max(x_f)

v_f=[i[1] for i in sol_d[-500:]]
v_max=np.max(v_f)

plt.figure(figsize=(8,6))
plt.plot(t_f,[i/x_max for i in x_f],label=r'$X$(t)')
plt.plot(t_f,[i/v_max for i in v_f],label=r'$\dot{X}$(t)')
plt.plot(t_f,fuerza_f,label=r'$\sin(\Theta$(t))')
plt.xlabel('Tiempo (s)',fontsize=label_size)
plt.yticks([],[])
plt.xlim([98,100])
plt.legend(fontsize=14,loc=3)
plt.grid('on')
plt.savefig('fases.png',dpi=300,bbox_inches='tight',facecolor='white')

