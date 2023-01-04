#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import math as mt
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root
from iapws import IAPWS97
import iapws
import fluids
import pandoc

### pi
pi=mt.atan(1.)*4.

### gravidade no Porto
g=fluids.gravity(latitude=41.15, H=50) 

#Propriedades do fluido
p_estipulado=0.200 #Pressão estipulada em MPa
water=IAPWS97(T=90+273.15, P=p_estipulado)
rho=water.Liquid.rho
mu=water.Liquid.mu
nu=water.Liquid.nu

#Criação de variavéis
mq_ptotal=[0,0,0,0,0]
mq_hi=[0,0,0,0,0]
mq_hf=[0,0,0,0,0]
mq_caudal=[0,0,0,0,0]
mq_caudalmass=[0,0,0,0,0]

########## DADOS ###############
mq_deltap=[60000,60000,60000,60000,60000] #Perda de carga das máquinas (Pa)
mq_ti=[90,90,90,90,90] #Temperatura inicial das máquinas (ºC)
mq_tf=[70,70,70,70,70] #Temperatura final das máquinas (ºC)

for x in range(0,5):          #Conversão de graus para kelvins
    mq_ti[x]=mq_ti[x]+273.15
    mq_tf[x]=mq_tf[x]+273.15

mq_putil=[38,188,263,563,450] #Potencia util das máquinas (kW) 
mq_eficiencia=[0.92,0.92,0.92,0.92,0.92] #Eficiencia das máquinas

#Cálculo Potência Util (kW)
for x in range(0,5):
    mq_ptotal[x]=mq_putil[x]/mq_eficiencia[x]

#Dimensões geométricas
fabrica_altura=5 #Altura da fábrica (m)

mq_comprimento=[3,4,5,7,6] #Comprimento das máquinas [m]
mq_largura=[3,3,4,4,4] #Largura das máquinas [m]
mq_altura=[2,2,2,3,2,2]
mq_espacamento=2 #Espaçamento entre máquinas [m]


# In[ ]:


####### CÁLCULO CAUDAL MÁSSICO ##########

for x in range (0,5): #Calculo da entalpia inicial e final em cada máquina
    water=IAPWS97(T=mq_ti[x], P=p_estipulado)
    mq_hi[x]=(float(water.h))
    water=IAPWS97(T=mq_tf[x], P=p_estipulado)
    mq_hf[x]=(float(water.h))
    

temp_caudaltotal=0
for x in range (0,5):
    mq_caudalmass[x]=mq_ptotal[x]/(mq_hi[x]-mq_hf[x]) #Cálculo do caudal Potencia total/(hi-hf)
    mq_caudal[x]=mq_caudalmass[x]/rho
    print("Máquina",x+1,":\n","Caudal mássico:",round(mq_caudalmass[x],3), "kg/s", "\n Caudal volumico:",round(mq_caudal[x]*3600,3),"m^3/h")
    temp_caudaltotal+=mq_caudalmass[x]

print("\nTotal:\nCaudal mássico:",round(temp_caudaltotal,3),"kg/s","\nCaudal volumico:",round(temp_caudaltotal/rho*3600,3),"m^3/h\n                ",round(temp_caudaltotal/rho*3600*1000,3),"l/h")


# In[ ]:


#Criação de variaveis
segmento_RC=[0,0,0,0,0,0]
segmentos_alimentacao=np.zeros((11,7),dtype=np.float64) #Matriz com Q,comprimento,Jlinha,k,Jlocalizada
segmentos_saida=np.zeros((11,7),dtype=np.float64) #Matriz com Q,comprimento,Jlinha,k,Jlocalizada
segmentos_jtotal=np.zeros((11,6),dtype=float) #Matriz de somatorio com JLinha, JLocalizada, JMaquinas, Diferença altura, e Jtotal
segmentos_designacoes=[0,"SS12345","SS345  ","SS12   ","SS1    ","SS2    ","SS3    ","SS4    ","SS5    ","SS34   ","SS45   "]

#Dados do tubo
eps=0.15/1000 #Rugosidade do tubo 

#Calculo de caudais por piso
piso1_caudalmasstotal=mq_caudalmass[2]+mq_caudalmass[3]+mq_caudalmass[4] 
piso2_caudalmasstotal=mq_caudalmass[0]+mq_caudalmass[1]
print("caudal total piso1:",piso1_caudalmasstotal/rho*3600,"m^3/h\ncaudalmass total piso2:",piso2_caudalmasstotal/rho*3600,"m^3/h")

######## CALCULO CAUDAIS E COMPRIMENTOS POR SEGMENTOS ########

#Calculos caudais por segmentos

segmento_RC[1]=mq_caudal[0]+mq_caudal[1]+mq_caudal[2]+mq_caudal[3]+mq_caudal[4] #RC

segmentos_alimentacao[1,1]=mq_caudal[0]+mq_caudal[1]+mq_caudal[2]+mq_caudal[3]+mq_caudal[4] #A12345
segmentos_alimentacao[2,1]=mq_caudal[2]+mq_caudal[3]+mq_caudal[4] #A345
segmentos_alimentacao[3,1]=mq_caudal[0]+mq_caudal[1] #A12
segmentos_alimentacao[4,1]=mq_caudal[0] #A1
segmentos_alimentacao[5,1]=mq_caudal[1] #A2
segmentos_alimentacao[6,1]=mq_caudal[2] #A3
segmentos_alimentacao[7,1]=mq_caudal[3] #A4
segmentos_alimentacao[8,1]=mq_caudal[4] #A5
segmentos_alimentacao[9,1]=mq_caudal[2]+mq_caudal[3] #A34
segmentos_alimentacao[10,1]=0 #A45 INEXISTENTE!!

#Definição de caudais por segmento

for x in range(1,11):
    if x in [1,2,3,4,5,6,7,8]:
        segmentos_saida[x,1]=segmentos_alimentacao[x,1] #Caudais iguais na alimentacao e na saida no mesmo segmento
    elif x==9: #S34 INEXISTENTE!!
        segmentos_saida[x,1]=0 
    elif x==10: #A45 INEXISTENTE!!
        segmentos_saida[x,1]=mq_caudal[3]+mq_caudal[4]
        
#Calculos comprimentos por segmentos

segmento_RC[2]=6.388 #RC

segmentos_alimentacao[1,2]=4.454 #A12345
segmentos_alimentacao[2,2]=3.747 #A345
segmentos_alimentacao[3,2]=2.459 #fabrica_altura+0.5 #A12
segmentos_alimentacao[4,2]=5.675 #mq_largura[1]/2+mq_espacamento+mq_largura[0]+0.5/2#A1
segmentos_alimentacao[5,2]=0.708 #A2
segmentos_alimentacao[6,2]=7.37 #mq_largura[2]/2+mq_espacamento+mq_largura[3]/2+0.5 #A3
segmentos_alimentacao[7,2]=1.549 #A4
segmentos_alimentacao[8,2]=7.643 #mq_largura[3]/2+mq_espacamento+mq_largura[4]/2+0.5 #A5
segmentos_alimentacao[9,2]=1.002 #A34
segmentos_alimentacao[10,2]=0 #A45 INEXISTENTE!!

segmentos_saida[1,2]=0.167 #S12345
segmentos_saida[2,2]=0.929 #S345
segmentos_saida[3,2]=6.406 #fabrica_altura+1 #S12
segmentos_saida[4,2]=2.692 #mq_largura[1]/2+mq_espacamento+mq_largura[0]/2+0.5#S1
segmentos_saida[5,2]=3.431 #S2
segmentos_saida[6,2]=2.306 #mq_largura[2]/2+mq_espacamento+mq_largura[3]/2+0.5 #S3
segmentos_saida[7,2]=1.246 #S4
segmentos_saida[8,2]=8.274 #mq_largura[3]/2+mq_espacamento+mq_largura[4]/2+0.5 #S5
segmentos_saida[9,2]=0 #S34 INEXISTENTE!!
segmentos_saida[10,2]=4.815 #S45


#############################################################
############# FUNCAO DE CALCULO DAS PERDAS ##################

def calculo_perdastotais(d1,d2,d3,d4,d5,d6,d7,d8):
    
    #RC,A12345,A345,A12,A1,A2,A3,A4,A5,A34,A45
    
    #Criação de variaveis
    
    d_segmentos=[d1,d1,d1,d6,d7,d6,d5,d3,d4,d1,d2]
    v_segmentos=[0,0,0,0,0,0,0,0,0,0,0]
    rey_segmentos=[0,0,0,0,0,0,0,0,0,0,0]
    
    comp_tubo=[segmento_RC[2],0,0,0,0,0,0,0,0]
 
    #Calculo comprimento do tubo
    
    for x in range (0,11):
        if d_segmentos[x]==d1:
            comp_tubo[1]+=segmentos_alimentacao[x,2]+segmentos_saida[x,2]
        elif d_segmentos[x]==d2:
            comp_tubo[2]+=segmentos_alimentacao[x,2]+segmentos_saida[x,2]
        elif d_segmentos[x]==d3:
            comp_tubo[3]+=segmentos_alimentacao[x,2]+segmentos_saida[x,2]
        elif d_segmentos[x]==d4:
            comp_tubo[4]+=segmentos_alimentacao[x,2]+segmentos_saida[x,2]
        elif d_segmentos[x]==d5:
            comp_tubo[5]+=segmentos_alimentacao[x,2]+segmentos_saida[x,2]
        elif d_segmentos[x]==d6:
            comp_tubo[6]+=segmentos_alimentacao[x,2]+segmentos_saida[x,2]
        elif d_segmentos[x]==d7:
            comp_tubo[7]+=segmentos_alimentacao[x,2]+segmentos_saida[x,2]
        elif d_segmentos[x]==d8:
            comp_tubo[8]+=segmentos_alimentacao[x,2]+segmentos_saida[x,2]
            
    #Cálculo Perdas em Linha por segmento
    
    for x in range (0,11):
        if x in [1,2,3,4,5,6,7,8,9]:
            v_segmentos[x]=segmentos_alimentacao[x,1]/(pi*d_segmentos[x]**2/4)
        elif x==10:
            v_segmentos[x]=segmentos_saida[x,1]/(pi*d_segmentos[x]**2/4)
        elif x==0:
            v_segmentos[x]=segmento_RC[1]/(pi*d_segmentos[x]**2/4)
            
        rey_segmentos[x]=fluids.Reynolds(V=v_segmentos[x], D=d_segmentos[x], rho=rho, mu=mu)
        f=fluids.friction.friction_factor(Re=rey_segmentos[x],eD=eps/d_segmentos[x],Method='Moody')
        if x==0:
            segmento_RC[3]=(f*segmento_RC[2]/d_segmentos[x])*((8*segmento_RC[1]**2)/(pi**2*d_segmentos[x]**4*g))
        else:
            segmentos_alimentacao[x,3]=(f*segmentos_alimentacao[x,2]/d_segmentos[x])*((8*segmentos_alimentacao[x,1]**2)/(pi**2*d_segmentos[x]**4*g))
            segmentos_saida[x,3]=(f*segmentos_saida[x,2]/d_segmentos[x])*(8*segmentos_saida[x,1]**2/(pi**2*d_segmentos[x]**4*g))

    #### Calculos Perdas Localizadas por segmento ####

    for x in range(1,11): #Fazer reset aos k
        segmentos_alimentacao[x,4]=0
        segmentos_saida[x,4]=0
    segmento_RC[4]=0

    #Segmento RC (4 curvas, 2 tees, 1 borboletas, 1 antiretorno, 1 saida, 1 entrada)
    segmento_RC[4]+=4*fluids.Hooper2K(Di=d_segmentos[0]/0.0254, Re=rey_segmentos[0], K1=800, Kinfty=0.25) #Cotovelo flangeado 90 graus
    segmento_RC[4]+=1
    segmento_RC[4]+=1
    segmento_RC[4]+=2*fluids.fittings.K_branch_diverging_Crane(d_segmentos[0], d_segmentos[0],segmento_RC[1],0, angle=90)
    segmento_RC[4]+=1*fluids.Hooper2K(Di=d_segmentos[0]/0.0254, Re=rey_segmentos[0], K1=800, Kinfty=0.25) # Válvula Borboleta (Butterfly)
    segmento_RC[4]+=1*fluids.Hooper2K(Di=d_segmentos[0]/0.0254, Re=rey_segmentos[0], K1=1500, Kinfty=1.5) #Valvula Anti-Retorno
    
    ### Circuito alimentacao ###
    #Segmento A12345 (1 curva,1 tee divergente)
    segmentos_alimentacao[1,4]+=0.78
    segmentos_alimentacao[1,4]+=1*fluids.Hooper2K(Di=d_segmentos[1]/0.0254, Re=rey_segmentos[1], K1=800, Kinfty=0.25) #Cotovelo flangeado 90 graus
    segmentos_alimentacao[1,4]+=fluids.fittings.K_branch_diverging_Crane(d_segmentos[1], d_segmentos[3],segmentos_alimentacao[1,1],segmentos_alimentacao[3,1], angle=90) #tee divergente
    
    #Segmento A345 (1 curva, 2 tee)
    segmentos_alimentacao[2,4]+=fluids.fittings.K_branch_diverging_Crane(d_segmentos[2], d_segmentos[8],segmentos_alimentacao[2,1],segmentos_alimentacao[8,1], angle=90)
    segmentos_alimentacao[2,4]+=fluids.fittings.K_branch_diverging_Crane(d_segmentos[2], d_segmentos[6],segmentos_alimentacao[2,1],segmentos_alimentacao[6,1], angle=90)
    segmentos_alimentacao[2,4]+=1*fluids.Hooper2K(Di=d_segmentos[2]/0.0254, Re=rey_segmentos[2], K1=800, Kinfty=0.25) #Cotovelo flangeado 90 graus
    
    #Segmento A12 (1 curva, 1 tee divergente)
    segmentos_alimentacao[3,4]+=1*fluids.Hooper2K(Di=d_segmentos[3]/0.0254, Re=rey_segmentos[3], K1=800, Kinfty=0.25) #Cotovelo flangeado
    segmentos_alimentacao[3,4]+=fluids.fittings.K_branch_diverging_Crane(d_segmentos[3], d_segmentos[4],segmentos_alimentacao[3,1],segmentos_alimentacao[4,1], angle=90) #tee divergente
    
    #Segmento A1 (1 curva, 1 entrada)
    segmentos_alimentacao[4,4]+=1 #Entrada
    segmentos_alimentacao[4,4]+=1*fluids.Hooper2K(Di=d_segmentos[4]/0.0254, Re=rey_segmentos[4], K1=800, Kinfty=0.25) #Cotovelo flangeado 90º
    
    #Segmento A2 (1 entrada)
    segmentos_alimentacao[5,4]+=1
    
    #Segmento A34 ()
    segmentos_alimentacao[9,4]+=1*fluids.Hooper2K(Di=d_segmentos[9]/0.0254, Re=rey_segmentos[9], K1=800, Kinfty=0.25) #Cotovelo flangeado 90º
    segmentos_alimentacao[9,4]+=fluids.fittings.K_branch_diverging_Crane(d_segmentos[9], d_segmentos[7],segmentos_alimentacao[9,1],segmentos_alimentacao[7,1], angle=90)
    
    #Segmento A3 (1 curva,1 entrada)
    segmentos_alimentacao[6,4]+=1
    segmentos_alimentacao[6,4]+=1*fluids.Hooper2K(Di=d_segmentos[6]/0.0254, Re=rey_segmentos[6], K1=800, Kinfty=0.25) #Cotovelo flangeado 90º
    
    #Segmento A4 (1 entrada)
    segmentos_alimentacao[7,4]+=1
    
    #Segmento A5 (1 curva, 1 entrada)
    segmentos_alimentacao[8,4]=1
    segmentos_alimentacao[8,4]=1*fluids.Hooper2K(Di=d_segmentos[8]/0.0254, Re=rey_segmentos[8], K1=800, Kinfty=0.25) #Cotovelo flangeado 90º
    
    ###Circuito saida ###
    
    #Segmento S12345 (1 tee convergente)
    segmentos_saida[1,4]+=fluids.fittings.K_branch_converging_Crane(d_segmentos[3], d_segmentos[2],segmentos_saida[3,1],segmentos_saida[2,1], angle=90)
    
    #Segmento S345 (1 curva, 1 tee convergente)
    segmentos_saida[2,4]+=1*fluids.Hooper2K(Di=d_segmentos[2]/0.0254, Re=rey_segmentos[2], K1=800, Kinfty=0.25) #Cotovelo flangeado 90º
    segmentos_saida[2,4]+=fluids.fittings.K_branch_converging_Crane(d_segmentos[2], d_segmentos[6],segmentos_saida[2,1],segmentos_saida[6,1], angle=90)
    
    #Segmento S45 (1 tee)
    segmentos_saida[10,4]+=fluids.fittings.K_branch_converging_Crane(d_segmentos[10], d_segmentos[7],segmentos_saida[10,1],segmentos_saida[7,1], angle=90)
    
    #Segmento S12 (2 curva, 1 tee divergente)
    segmentos_saida[3,4]+=2*fluids.Hooper2K(Di=d_segmentos[3]/0.0254, Re=rey_segmentos[3], K1=800, Kinfty=0.25) #Cotovelo flangeado
    segmentos_saida[3,4]+=fluids.fittings.K_branch_converging_Crane(d_segmentos[3], d_segmentos[4],segmentos_saida[3,1],segmentos_saida[4,1], angle=90)
    
    #Segmento S1 (1 curva, 1 butterfly, 1 check, 1 saida)
    segmentos_saida[4,4]+=1 #Saida
    segmentos_saida[4,4]+=1*fluids.Hooper2K(Di=d_segmentos[4]/0.0254, Re=rey_segmentos[4], K1=1500, Kinfty=1.5) # Valvula Antiretorno (Check Swing)
    segmentos_saida[4,4]+=1*fluids.Hooper2K(Di=d_segmentos[4]/0.0254, Re=rey_segmentos[4], K1=800, Kinfty=0.25) #Cotovelo flangeado 90º
    segmentos_saida[4,4]+=1*fluids.Hooper2K(Di=d_segmentos[4]/0.0254, Re=rey_segmentos[4], K1=800, Kinfty=0.25) # Válvula Borboleta (Butterfly)
    
    #Segmento S2 (1 curva, 1 butterfly, 1 check, 1 saida)
    segmentos_saida[5,4]+=1
    segmentos_saida[5,4]+=1*fluids.Hooper2K(Di=d_segmentos[5]/0.0254, Re=rey_segmentos[5], K1=1500, Kinfty=1.5) # Valvula Antiretorno (Check Swing)
    segmentos_saida[5,4]+-1*fluids.Hooper2K(Di=d_segmentos[5]/0.0254, Re=rey_segmentos[5], K1=800, Kinfty=0.25) #Cotovelo flangeado 90º
    segmentos_saida[5,4]+-1*fluids.Hooper2K(Di=d_segmentos[5]/0.0254, Re=rey_segmentos[5], K1=800, Kinfty=0.25) #Válvula Borboleta (Butterfly)
    
    #Segmento S3 (3 curvas, 1 butterfly, 1 check, 1 saida)
    segmentos_saida[6,4]+=1
    segmentos_saida[6,4]+=1*fluids.Hooper2K(Di=d_segmentos[6]/0.0254, Re=rey_segmentos[6], K1=1500, Kinfty=1.5) # Valvula Antiretorno (Check Swing)
    segmentos_saida[6,4]+=3*fluids.Hooper2K(Di=d_segmentos[6]/0.0254, Re=rey_segmentos[6], K1=800, Kinfty=0.25) #Cotovelo flangeado 90º
    segmentos_saida[6,4]+=1*fluids.Hooper2K(Di=d_segmentos[6]/0.0254, Re=rey_segmentos[6], K1=800, Kinfty=0.25) #Valvula Borboleta (Butterfly)
    
    #Segmento S4 (2 curvas, 1 butterfly, 1 check, 1 saida)
    segmentos_saida[7,4]+=1
    segmentos_saida[7,4]+=2*fluids.Hooper2K(Di=d_segmentos[7]/0.0254, Re=rey_segmentos[7], K1=800, Kinfty=0.25) #Cotovelo
    segmentos_saida[7,4]+=1*fluids.Hooper2K(Di=d_segmentos[7]/0.0254, Re=rey_segmentos[7], K1=1500, Kinfty=1.5) # Valvula Antiretorno (Check Swing)
    segmentos_saida[7,4]+=1*fluids.Hooper2K(Di=d_segmentos[7]/0.0254, Re=rey_segmentos[7], K1=800, Kinfty=0.25) #Valvula Borboleta (Butterfly)
    
    #Segmento S5 (3 curvas, 1 butterfly, 1 check, 1 saida)
    segmentos_saida[8,4]+=1
    segmentos_saida[8,4]+=1*fluids.Hooper2K(Di=d_segmentos[4]/0.0254, Re=rey_segmentos[4], K1=1500, Kinfty=1.5) # Valvula Antiretorno (Check Swing)
    segmentos_saida[8,4]=3*fluids.Hooper2K(Di=d_segmentos[8]/0.0254, Re=rey_segmentos[8], K1=800, Kinfty=0.25) #Cotovelo flangeado 90º
    segmentos_saida[8,4]=1*fluids.Hooper2K(Di=d_segmentos[8]/0.0254, Re=rey_segmentos[8], K1=800, Kinfty=0.25) #Válvula borboleta (Butterfly)
    
    for x in range (1,11):
        segmentos_alimentacao[x,5]=(segmentos_alimentacao[x,4])*(8*segmentos_alimentacao[x,1]**2/(pi**2*d_segmentos[x]**4*g))
        segmentos_saida[x,5]=(segmentos_saida[x,4])*(8*segmentos_saida[x,1]**2/(pi**2*d_segmentos[x]**4*g))
    
    segmento_RC[5]=(segmento_RC[4])*(8*segmento_RC[1]**2/(pi**2*d_segmentos[0]**4*g))
    
    
    ################# CALCULO E LÓGICA ############################
    
    for x in range(1,11):
        
        segmentos_alimentacao[x,6]=segmentos_alimentacao[x,5]+segmentos_alimentacao[x,3]
        segmentos_saida[x,6]=segmentos_saida[x,5]+segmentos_saida[x,3]
        segmentos_jtotal[x,1]=segmentos_alimentacao[x,3]+segmentos_saida[x,3]
        segmentos_jtotal[x,2]=segmentos_alimentacao[x,5]+segmentos_saida[x,5]
        
        if x in [4,5,6,7,8]:
            segmentos_jtotal[x,3]=mq_deltap[1]/(rho*g)
        elif x==3:
            segmentos_jtotal[x,4]=fabrica_altura-fabrica_altura
        segmentos_jtotal[x,5]=segmentos_jtotal[x,1]+segmentos_jtotal[x,2]+segmentos_jtotal[x,3]+segmentos_jtotal[x,4]
        
    if ((segmentos_jtotal[8,5]+segmentos_jtotal[10,5])>(segmentos_jtotal[7,5]+segmentos_jtotal[10,5])+segmentos_jtotal[9,5]) and ((segmentos_jtotal[8,5]+segmentos_jtotal[10,5])>(segmentos_jtotal[6,5]+segmentos_jtotal[9,5])): #Comparaçao entre 3 segmentos
        piso1_jtotal=segmentos_jtotal[1]+segmentos_jtotal[2]+segmentos_jtotal[8]+segmentos_jtotal[10]
        print("Perdas maiores pelos segmentos da Máquina 5 no Piso 1")
    elif (segmentos_jtotal[7,5]+segmentos_jtotal[10,5])>segmentos_jtotal[6,5]: #Perdas maq A4 > A3
        piso1_jtotal=segmentos_jtotal[1]+segmentos_jtotal[2]+segmentos_jtotal[9]+segmentos_jtotal[7]+segmentos_jtotal[10]
        print("Perdas maiores pelos segmentos da Máquina 4 no Piso 1")
    else:
        piso1_jtotal=segmentos_jtotal[1]+segmentos_jtotal[2]+segmentos_jtotal[9]+segmentos_jtotal[6]
        print("Perdas maiores pelos segmentos da Máquina 3 no Piso 1")
        
    
    if segmentos_jtotal[4,5]>segmentos_jtotal[5,5]: #Comparação de perdas maiores entre Segmento A1+S1 e Segmento A2+S2
        piso2_jtotal=segmentos_jtotal[1]+segmentos_jtotal[3]+segmentos_jtotal[4]
        print("Perdas maiores pelos segmentos da Máquina 1 no Piso 2")
    else:
        piso2_jtotal=segmentos_jtotal[1]+segmentos_jtotal[3]+segmentos_jtotal[5]
        print("Perdas maiores pelos segmentos da Máquina 2 no Piso 2\n")
    
    return (d_segmentos,v_segmentos,piso1_jtotal,piso2_jtotal,segmento_RC,segmentos_jtotal,segmentos_alimentacao,segmentos_saida,comp_tubo,rey_segmentos)
   


# In[ ]:


#Cálculo das perdas de carga para um diâmetro específico
d1=94.4/1000 #Em metros
d2=82.5/1000  #Em metros
d3=69.7/1000 #Em metros
d4=54.5/1000 #Em metros
d5=42.5/1000
d6=37.2/1000
d7=17.2/1000
d8=0/1000

d_segmentos,v_segmentos,piso1_jtotal,piso2_jtotal,segmento_RC,segmentos_jtotal,segmentos_alimentacao,segmentos_saida,comp_tubo,rey_segmentos=calculo_perdastotais(d1,d2,d3,d4,d5,d6,d7,d8)
print("\nPara o diâmetro:")
for x in range (1,9):
    print("Comprimento d",x,":",round(comp_tubo[x],3),"m")
    
if piso1_jtotal[5] > piso2_jtotal[5]:
    print("   \nPerdas totais:",round(piso1_jtotal[5]+segmento_RC[3]+segmento_RC[5],3),"m\n")
else:
    print("   \nPerdas totais:",round(piso2_jtotal[5]+segmento_RC[3]+segmento_RC[5],3),"m\n")
print("Perdas Piso 1:",round(piso1_jtotal[5],3),"m   Piso 2:",round(piso2_jtotal[5],3),"m   RC:",round(segmento_RC[3]+segmento_RC[5],3),"m")
print("Velocidade por segmentos Piso 1:","A12345:",round(v_segmentos[1],3),"m/s  A345",round(v_segmentos[2],3),"m/s  A3",round(v_segmentos[6],3),"m/s","m/s  A4",round(v_segmentos[7],3),"m/s","m/s  A5",round(v_segmentos[8],3),"m/s")
print("Velocidade por segmentos Piso 2:","A12:",round(v_segmentos[3],3),"m/s  A1",round(v_segmentos[4],3),"m/s","m/s  A2",round(v_segmentos[5],3),"m/s")

print("\nCircuito Alimentacao\n    -   |    Q    |    L    |  J Linha |    k    | J localizada  | J total")
for x in range (1,11):
    print(segmentos_designacoes[x],"",round(segmentos_alimentacao[x,1],6),"   ",round(segmentos_alimentacao[x,2],2),"   ",round(segmentos_alimentacao[x,3],4),"    ",round(segmentos_alimentacao[x,4],2),"    ",round(segmentos_alimentacao[x,5],3),"        ",round(segmentos_alimentacao[x,6],3))
print("\nCircuito Saida   \n    -   |   Q    |   L  |  J Linha |    k    | J localizada  | J total")
for x in range (1,11):
    print(segmentos_designacoes[x],"",round(segmentos_saida[x,1],4),"   ",round(segmentos_saida[x,2],2),"   ",round(segmentos_saida[x,3],4),"    ",round(segmentos_saida[x,4],2),"    ",round(segmentos_saida[x,5],3),"        ",round(segmentos_saida[x,6],3))
print("\nCircuito RC   \n    -    |   Q  |   L  |  J Linha |    k    | J localizada  | J total")
print("RC       ",round(segmento_RC[1],4)," ",round(segmento_RC[2],4),"  ",round(segmento_RC[3],4),"  ",round(segmento_RC[4],4),"   ",round(segmento_RC[5],4),"    ",round(segmento_RC[3]+segmento_RC[5],4))
#Calculo velocidade minima antiretorno
water=IAPWS97(T=70+273.15, P=p_estipulado)
rho70=water.Liquid.rho


# In[ ]:


import ht
from numpy import log as ln
import sympy as sym

#Calculo térmico
stefan_boltzmann=5.67*10**-8
emiss=0

p_estipulado=0.200 #Pressão estipulada em MPa
t_atm=20
h2_arbitado=100

esp=0/1000 #espessura do isolante

k2=56.7
k3=0.04 #W/(m.K)

n_div=20

#################################################################
#################  DEFINIÇÃO DE VARIÁVEIS  ######################

segmentos_alimentacao_termico=np.zeros((11,10),dtype=np.float64) #Matriz com L,d1,d2,e,rey,t_inicial,t_final,q
segmentos_saida_termico=np.zeros((11,10),dtype=np.float64) #Matriz com  L,d1,d2,e,rey,t_inicial,t_final,q
segmento_RC_termico=[0,0,0,0,0,0,0,0,0,0]
d2_segmentos=[d1+3.6*2/1000,d1+3.6*2/1000,d1+3.6*2/1000,d6+2.6*2/1000,d7+2*2/1000,d6+2.6*2/1000,d5+2.9*2/1000,d3+3.2*2/1000,d4+2.9*2/1000,d1+3.6*2/1000,d2+3.2*2/1000]
#print(d2_segmentos)
#Transposição de dados de variaveis do calculo de perdas para a térmica
water=IAPWS97(T=90+273.15, P=p_estipulado)

segmento_RC_termico[1]=segmento_RC[1]
segmento_RC_termico[2]=d_segmentos[0]
segmento_RC_termico[3]=d2_segmentos[0]
segmento_RC_termico[5]=rey_segmentos[0]
segmento_RC_termico[9]=segmento_RC[1]*water.Liquid.rho

for x in range(1,11):
    segmentos_alimentacao_termico[x,1]=segmentos_alimentacao[x,2]
    segmentos_saida_termico[x,1]=segmentos_saida[x,2]
    
    segmentos_alimentacao_termico[x,2]=d_segmentos[x]
    segmentos_alimentacao_termico[x,3]=d2_segmentos[x]
    segmentos_saida_termico[x,2]=d_segmentos[x]
    segmentos_saida_termico[x,3]=d2_segmentos[x]
    
    segmentos_alimentacao_termico[x,5]=rey_segmentos[x]
    segmentos_saida_termico[x,5]=rey_segmentos[x]
    
    segmentos_alimentacao_termico[x,9]=segmentos_alimentacao[x,1]*water.Liquid.rho
    segmentos_saida_termico[x,9]=segmentos_saida[x,1]*water.Liquid.rho

    segmentos_alimentacao_termico[x,4]=esp
    segmentos_saida_termico[x,4]=esp
    
#################################################################
###########################  FUNCÕES  ###########################

def call_ciclos(matriz,linha):
    return calculo_principal(matriz[linha,1],matriz[linha,2],matriz[linha,3],matriz[linha,4],matriz[linha,5],matriz[linha,6],matriz[linha,9])

def calculo_principal(comprimento,d1,d2,esp,rey,t_inicial,caudalmass):
    
    #Definição de variaveis
    global n_div
    x_varr=0
    q_total=0
    deltax=comprimento/n_div
    t=t_inicial
    
    for x in range(1,n_div+1):
        t_saida,q_inst=fluxo_instantaneo(t,deltax,d1,d2,esp,rey,caudalmass)
        q_total+=q_inst
        t=t_saida
    return(t,q_total)

def fluxo_instantaneo(t_inicial,l,d1,d2,esp,rey,m):
    
    global t_atm,k2,k3
    
    water=IAPWS97(T=t_inicial+273.15, P=p_estipulado)
    prandt=water.Liquid.Prandt
    cp=water.Liquid.cp*1000
    k1=water.Liquid.k  
    nu=ht.conv_internal.turbulent_Dittus_Boelter(rey,prandt, heating=False, revised=True) #Equação Dittus-Boelter, arrefecimento e com m=0.023
    h1=nu*k1/d1 #Calculo h interno
    r1=1/(h1*pi*d1*l) #Calculo resistencia da convecção interna
    r2=(ln(d2/d1))/(2*pi*k2*l) #Calculo da condução do tubo
    d3=d2+esp*2 #Diametro exterior do isolante
    r3=(ln(d3/d2))/(2*pi*k3*l) #Calculo da condução no isolante

    n_iter=1
    erro=1
    h2_arbitado=100
    
    while erro > 0.001 or erro <-0.001:
        #print("\nIteração",n_iter,"\nerro:",erro,"\nh2:",h2_arbitado,"Nuext:",nu)
        r_paralelo=1/(h2_arbitado*pi*d3*l)
        q_inst=(t_inicial-t_atm)/(r1+r2+r3+r_paralelo)

        t1=t_inicial-q_inst*(r1)
        t2=t_inicial-q_inst*(r1+r2)
        t3=t_inicial-q_inst*(r1+r2+r3)

        #print("t1:",round(t1,3),"t2:",round(t2,3),"t3:",round(t3,3),"\nFluxo de calor:",round(q_inst,3),"W")

        air=iapws.humidAir.Air(T=((t_atm+t3)/2+273.15),P=0.101325)
        visc=air.nu #viscosidade cinemática
        prandt=air.Prandt #numero de Prandt
        k4=air.k #k ar
        b=1/((t_atm+t3)/2+273.15) #Coeficiente de expansão- gás ideal

        gr=(g*b*(t3-t_atm)*d3**2)/visc**2
        nu=ht.conv_free_immersed.Nu_horizontal_cylinder_Churchill_Chu(prandt,gr)
        h2_conv=nu*k4/d3
        h2_rad=emiss*stefan_boltzmann*((t3+273.15)**2+(t_atm+273.15)**2)*(t_atm+273.15+t3+273.15)
        h2=h2_conv+h2_rad

        erro=1-h2_arbitado/h2
        n_iter=n_iter+1
        h2_arbitado=h2
        #print(h2)

    x=sym.symbols('x')
    eq1 = sym.Eq(q_inst,m*cp*(t_inicial-x))
    t_saida=float(sym.solve((eq1),(x))[0])
    #print("\nTemperatura saida:",round(t_saida,8),"ºC")
    return(t_saida,q_inst)

def mistura_fluidos(matriz,l1,l2):

    m1,m2=matriz[l1,9],matriz[l2,9]
    t1,t2=matriz[l1,7],matriz[l2,7]
    
    water=IAPWS97(T=t1+273.15, P=p_estipulado)
    cp1=water.Liquid.cp*1000
    water=IAPWS97(T=t2+273.15, P=p_estipulado)
    cp2=water.Liquid.cp*1000
    
    x=sym.symbols('x')
    eq1 = sym.Eq(m1*cp1*(x-t1),-m2*cp2*(x-t2))
    mx=sym.solve(eq1)
    return sym.solve((eq1),(x))[0]

#################################################################
#######################  CÁLCULO E LÓGICA  ######################

#Circuito de entrada

t_inicial_alimentacao=90

segmentos_alimentacao_termico[1,6]=t_inicial_alimentacao
segmentos_alimentacao_termico[1,7],segmentos_alimentacao_termico[1,8]=call_ciclos(segmentos_alimentacao_termico,1)
segmentos_alimentacao_termico[2,6]=segmentos_alimentacao_termico[1,7] #Temp. entrada A345 = Temp. saida A12345
segmentos_alimentacao_termico[3,6]=segmentos_alimentacao_termico[1,7] #Temp. entrada A12 = Temp. saida A12345
segmentos_alimentacao_termico[2,7],segmentos_alimentacao_termico[2,8]=call_ciclos(segmentos_alimentacao_termico,2) #Calculo A345
segmentos_alimentacao_termico[3,7],segmentos_alimentacao_termico[3,8]=call_ciclos(segmentos_alimentacao_termico,3) #Calculo A12
segmentos_alimentacao_termico[4,6]=segmentos_alimentacao_termico[3,7] #Temp. entrada A1 = Temp. saida A12
segmentos_alimentacao_termico[5,6]=segmentos_alimentacao_termico[3,7] #Temp. entrada A2 = Temp. saida A12
segmentos_alimentacao_termico[4,7],segmentos_alimentacao_termico[4,8]=call_ciclos(segmentos_alimentacao_termico,4) #Calculo A1
segmentos_alimentacao_termico[5,7],segmentos_alimentacao_termico[5,8]=call_ciclos(segmentos_alimentacao_termico,5) #Calculo A2
segmentos_alimentacao_termico[9,6]=segmentos_alimentacao_termico[2,7] #Temp. entrada A34 = Temp. saida A345
segmentos_alimentacao_termico[8,6]=segmentos_alimentacao_termico[2,7] #Temp. entrada A5 = Temp. saida A345
segmentos_alimentacao_termico[9,7],segmentos_alimentacao_termico[9,8]=call_ciclos(segmentos_alimentacao_termico,9) #Calculo A34
segmentos_alimentacao_termico[8,7],segmentos_alimentacao_termico[8,8]=call_ciclos(segmentos_alimentacao_termico,8) #Calculo A5
segmentos_alimentacao_termico[6,6]=segmentos_alimentacao_termico[9,7] #Temp. entrada A3 = Temp. saida A34
segmentos_alimentacao_termico[7,6]=segmentos_alimentacao_termico[9,7] #Temp. entrada A4 = Temp. saida A34
segmentos_alimentacao_termico[6,7],segmentos_alimentacao_termico[6,8]=call_ciclos(segmentos_alimentacao_termico,6) #Calculo A3
segmentos_alimentacao_termico[7,7],segmentos_alimentacao_termico[7,8]=call_ciclos(segmentos_alimentacao_termico,7) #Calculo A4

#Circuito de Saida

t_inicial_saida=70
segmentos_saida_termico[4,6]=t_inicial_saida #Temp. Inicial A1
segmentos_saida_termico[5,6]=t_inicial_saida #Temp. Inicial A2
segmentos_saida_termico[4,7],segmentos_saida_termico[4,8]=call_ciclos(segmentos_saida_termico,4) #Calculo A1
segmentos_saida_termico[5,7],segmentos_saida_termico[5,8]=call_ciclos(segmentos_saida_termico,5) #Calculo A2
segmentos_saida_termico[3,6]=mistura_fluidos(segmentos_saida_termico,4,5) #Mistura entre A1 e A2 - Temp. Inicial A12
segmentos_saida_termico[3,7],segmentos_saida_termico[3,8]=call_ciclos(segmentos_saida_termico,3) #Calculo A12
segmentos_saida_termico[6,6]=t_inicial_saida #Temp. Inicial A3
segmentos_saida_termico[7,6]=t_inicial_saida #Temp. Inicial A4
segmentos_saida_termico[8,6]=t_inicial_saida #Temp. Inicial A5
segmentos_saida_termico[6,7],segmentos_saida_termico[6,8]=call_ciclos(segmentos_saida_termico,6) #Calculo A3
segmentos_saida_termico[7,7],segmentos_saida_termico[7,8]=call_ciclos(segmentos_saida_termico,7) #Calculo A4
segmentos_saida_termico[8,7],segmentos_saida_termico[8,8]=call_ciclos(segmentos_saida_termico,8) #Calculo A5
segmentos_saida_termico[10,6]=mistura_fluidos(segmentos_saida_termico,7,8) #Mistura entre A5 e A4 - Temp. Inicial A45
segmentos_saida_termico[10,7],segmentos_saida_termico[10,8]=call_ciclos(segmentos_saida_termico,10) #Calculo A45
segmentos_saida_termico[2,6]=mistura_fluidos(segmentos_saida_termico,9,6) #Mistura entre A45 e A3 - Temp. Inicial A345
segmentos_saida_termico[2,7],segmentos_saida_termico[2,8]=call_ciclos(segmentos_saida_termico,2) #Calculo A345
segmentos_saida_termico[1,6]=mistura_fluidos(segmentos_saida_termico,3,2) #Mistura entre A345 e A12 - Temp. Inicial A12345
segmentos_saida_termico[1,7],segmentos_saida_termico[1,8]=call_ciclos(segmentos_saida_termico,1) #Calculo A12345

q_total_alimentacao=0
for x in range(1,11):
    q_total_alimentacao+=segmentos_alimentacao_termico[x,8]+segmentos_saida_termico[x,8]


# In[ ]:


print("Perdas Caloríficas Circuito Alimentação:",round(q_total_alimentacao,3), "W")

print("\nCircuito Alimentacao\n    -   |   L   |   d1   |   d2   |  esp  |    rey     |   t_inicial  |    t_final    |    q    |   caudal mass")
for x in range (1,11):
    print(segmentos_designacoes[x]," ",round(segmentos_alimentacao_termico[x,1],6),"  ",round(segmentos_alimentacao_termico[x,2],2),"  ",round(segmentos_alimentacao_termico[x,3],4),"  ",round(segmentos_alimentacao_termico[x,4],2)," ",round(segmentos_alimentacao_termico[x,5],3),"     ",round(segmentos_alimentacao_termico[x,6],5),"      ",round(segmentos_alimentacao_termico[x,7],5),"     ",round(segmentos_alimentacao_termico[x,8],3),"     ",round(segmentos_alimentacao_termico[x,9],3))
    
print("\nCircuito Saida   \n    -   |   L   |   d1   |   d2   |  esp  |    rey     |   t_inicial  |    t_final    |    q    |   caudal mass")
for x in range (1,11):
    print(segmentos_designacoes[x]," ",round(segmentos_saida_termico[x,1],6),"  ",round(segmentos_saida_termico[x,2],2),"  ",round(segmentos_saida_termico[x,3],4),"  ",round(segmentos_saida_termico[x,4],4)," ",round(segmentos_saida_termico[x,5],3),"     ",round(segmentos_saida_termico[x,6],5),"      ",round(segmentos_saida_termico[x,7],5),"     ",round(segmentos_saida_termico[x,8],3),"     ",round(segmentos_saida_termico[x,9],3))


# In[ ]:


######################################################################
########## OTIMIZACAO DE ESPESSURA DE ISOLAMENTO #####################

from scipy import optimize

segmentos_alimentacao_termico_opt=np.zeros((11,3),dtype=np.float64) #Matriz com esp_opt, custo_total
segmentos_saida_termico_opt=np.zeros((11,3),dtype=np.float64) #Matriz com esp_opt, custo_total

ins_cost=0

esp_variavel=0.05

fuel_price=0.0023 #[€/kWh]
op_time=365*24
eff=0.85 #heater efficiency

def variacao_esp(esp_variavel):
    global esp, segmentos_alimentacao_termico
    esp=esp_variavel
    segmentos_alimentacao_termico[linha,4]=esp
    segmentos_saida_termico[linha,4]=esp
    
    matriz[linha,7],matriz[linha,8]=call_ciclos(matriz,linha)
    
    fuel_cost=matriz[linha,8]*0.001*op_time/eff*fuel_price
    
    #Escolha de equações de custo do Isolamento
    if round(matriz[linha,3],4)==0.1016:
        ins_cost=(2*(1878.8*esp-7.42)/25)*matriz[linha,1]
    elif round(matriz[linha,3],4)==0.0889:
        ins_cost=(2*(1763.6*esp - 4.1355)/25)*matriz[linha,1]
    elif round(matriz[linha,3],4)==0.0761:
        ins_cost=(2*(1445.2*esp - 1.4767)/25)*matriz[linha,1]
    elif round(matriz[linha,3],4)==0.0603:
        ins_cost=(2*(1483.8*esp - 4.7392)/25)*matriz[linha,1]
    elif round(matriz[linha,3],4)==0.0483:
        ins_cost=(2*(1163.7*esp + 1.4054)/25)*matriz[linha,1]    
    elif round(matriz[linha,3],4)==0.0424:
        ins_cost=(2*(1195.4*esp - 3.5999)/25)*matriz[linha,1]
    elif round(matriz[linha,3],4)==0.0212:
        ins_cost=(2*(1141.5*esp - 7.0088)/25)*matriz[linha,1]
    else:
        print("Erro! - Sem Equação de Custo de Isolamento")
    return(ins_cost+fuel_cost)

def grafico(matriz,linha):
    x=[]
    y=[]
    esp=0/1000 #Iniciação do range (em m)
    while esp<(10/1000):
        linha=1
        matriz=segmentos_alimentacao_termico
        custo=variacao_esp(esp)
        plt.figure(1)
        plt.plot([esp*1000],[custo],'ro',label="piso1")
        esp=esp+(1/1000)
    plt.show()
    return()

def otimizacao(matriz_call,linha_call):
    global matriz, linha
    matriz=matriz_call
    linha=linha_call
    res = optimize.minimize(variacao_esp,esp, method='Nelder-Mead',tol=1e-2)
    return float(res.x), float(res.fun)

#################################################################
#######################  CÁLCULO E LÓGICA  ######################

#Circuito de entrada

segmentos_alimentacao_termico_opt[1,1],segmentos_alimentacao_termico_opt[1,2]=otimizacao(segmentos_alimentacao_termico,1)
segmentos_alimentacao_termico_opt[2,1],segmentos_alimentacao_termico_opt[2,2]=otimizacao(segmentos_alimentacao_termico,2)
segmentos_alimentacao_termico_opt[3,1],segmentos_alimentacao_termico_opt[3,2]=otimizacao(segmentos_alimentacao_termico,3)
segmentos_alimentacao_termico_opt[4,1],segmentos_alimentacao_termico_opt[4,2]=otimizacao(segmentos_alimentacao_termico,4)
segmentos_alimentacao_termico_opt[5,1],segmentos_alimentacao_termico_opt[5,2]=otimizacao(segmentos_alimentacao_termico,5)
segmentos_alimentacao_termico_opt[6,1],segmentos_alimentacao_termico_opt[6,2]=otimizacao(segmentos_alimentacao_termico,6)
segmentos_alimentacao_termico_opt[7,1],segmentos_alimentacao_termico_opt[7,2]=otimizacao(segmentos_alimentacao_termico,7)
segmentos_alimentacao_termico_opt[8,1],segmentos_alimentacao_termico_opt[8,2]=otimizacao(segmentos_alimentacao_termico,8)
segmentos_alimentacao_termico_opt[9,1],segmentos_alimentacao_termico_opt[9,2]=otimizacao(segmentos_alimentacao_termico,9)


#Circuito de Saida
segmentos_saida_termico_opt[1,1],segmentos_saida_termico_opt[1,2]=otimizacao(segmentos_saida_termico,1)
segmentos_saida_termico_opt[2,1],segmentos_saida_termico_opt[2,2]=otimizacao(segmentos_saida_termico,2)
segmentos_saida_termico_opt[3,1],segmentos_saida_termico_opt[3,2]=otimizacao(segmentos_saida_termico,3)
segmentos_saida_termico_opt[4,1],segmentos_saida_termico_opt[4,2]=otimizacao(segmentos_saida_termico,4)
segmentos_saida_termico_opt[5,1],segmentos_saida_termico_opt[5,2]=otimizacao(segmentos_saida_termico,5)
segmentos_saida_termico_opt[6,1],segmentos_saida_termico_opt[6,2]=otimizacao(segmentos_saida_termico,6)
segmentos_saida_termico_opt[7,1],segmentos_saida_termico_opt[7,2]=otimizacao(segmentos_saida_termico,7)
segmentos_saida_termico_opt[8,1],segmentos_saida_termico_opt[8,2]=otimizacao(segmentos_saida_termico,8)
segmentos_saida_termico_opt[10,1],segmentos_saida_termico_opt[10,2]=otimizacao(segmentos_saida_termico,10)


# In[ ]:


print("\nCircuito Alimentacao\n    -   | Espessura otimizada (mm) |   Custo Total  (€/ano)   |" )
for x in range (1,11):
    print(segmentos_designacoes[x],"          ",round(segmentos_alimentacao_termico_opt[x,1]*1000,2),"                ",round(segmentos_alimentacao_termico_opt[x,2],2))

print("\nCircuito Saida\n    -   | Espessura otimizada (mm) |   Custo Total  (€/ano)   |" )
for x in range (1,11):
    print(segmentos_designacoes[x],"          ",round(segmentos_saida_termico_opt[x,1]*1000,2),"                ",round(segmentos_saida_termico_opt[x,2],2))

