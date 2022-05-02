import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares,leastsq
import os
import argparse

def angle(x:complex)->float:
    return np.arctan(np.imag(x)/np.real(x))#x.imag/x.real)

def Maquina(Vsim,s,Wsincrono,Parametros):
    p   =   Parametros[0]
    Xm  =   Parametros[1]
    Rfe =   Parametros[2]
    R1  =   Parametros[3]
    R2  =   Parametros[4]
    X1  =   Parametros[5]
    X2  =   Parametros[6]

    j = complex(0,1)
    Zmag  = Rfe*j*Xm/(Rfe+j*Xm)
    Z2    = R2/s+j*X2
    Z1    = R1+j*X1
    Ztot  = Z1+Z2*Zmag/(Z2+Zmag)
    I1    = Vsim/Ztot #corrente nominal = Vnom/Ztot
    E     = Vsim-Z1*I1
    I2    = E/Z2
    modI2 = abs(I2)  

    Telesim = 3*(R2/s*modI2**2)/Wsincrono
    Iasim   = I1

    return Telesim, Iasim

def teste_maquina():
    T,I = Maquina(110,2,123,[1,2,3,4,5,6,7])
    print(f"T = {T}")
    print(f"I = {I}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', help='increase output verbosity', action='store_true')
    parser.add_argument('-s', '--save', help='save graphics in img folder', action='store_true')
    args = parser.parse_args()

    file_path = "MS_04_04_2022.csv"
    base_path = os.path.relpath(os.path.join(os.path.dirname(__file__),".."))
    img_path = os.path.join(base_path,"img")
    data_path = os.path.join(base_path,"data")
    data_file = os.path.join(data_path,file_path)


    print('+======= Experimento 4: ')
    # df = pd.read_csv(data_file, decimal=",")
    

    # entre aqui os parâmetros da sua máquina
    Vfase = 220/np.sqrt(3) # esta é a tensão de fase do modelo.
    frequencia = 60#
    p = 2 # numero de par de polos
    We = 2*np.pi*frequencia # frequência angular das grandezas elétricas rad/s
    Ws = We/p # esta é a rotação síncrona da máquina em rad/s
    RPM=frequencia*60/p  # esta é a rotação síncrona da máquina, em RPM

    # Parâmetros a serem introduzidos
    Xm=57.322
    Rf=374.22
    R1=2.7
    R2=3.9
    X1=2.73 # Aula
    X2=X1
    Par=[p, Xm, Rf, R1, R2, X1, X2]


    print("Escorregamento")
    # numero de pontos de escorregamento
    Ns=1000

    escorregamento = []
    # escorregamentos que serao calculados
    escormax=1
    for x in range(1,Ns+1):
        escr = (Ns-x)/(Ns-1)*escormax
        if escr == 0:
            escorregamento.append(1e-12) # não calcula para s=0
        else:
            escorregamento.append(escr)

    print("Curva")
    PerdasMecanicas = []
    rotacao = []
    Tele = []
    Ia = []
    Pmec = []
    Pentrada = []
    Rendimento = []
    # calculo da curva a partir do modelo
    for x in range(0,Ns): # aqui não é necessario fazer range[1,Ns+1], pois não será realizada operações com x, apenas será utilizado para indexação no array
        PerdasMecanicas.append(0)
        Wr = (1-escorregamento[x])*Ws
        rotacao.append(RPM*(1-escorregamento[x]))
        tele,ia = Maquina(Vfase,escorregamento[x],Ws,Par)    
        Tele.append(tele)    
        Ia.append(ia)    

        Pmec.append(Tele[x]*Wr)

        val = Vfase*Ia[x]
        Pentrada.append(3*val.real)
        Rendimento.append((Pmec[x]-PerdasMecanicas[x])/Pentrada[x])
    
    
    ModIa=np.abs(Ia)
    cosfi=np.cos(angle(Ia)) # a fase da tensão é a referência

    # # graficos em funcao da velocidade angular(pu)
    plt.figure()
    plt.subplot(2,1,1)
    plt.plot(rotacao,ModIa)
    plt.title("Ia(A) x rotação(rpm)")
    plt.subplot(2,1,2)
    plt.plot(rotacao,cosfi)
    plt.title("cos(fi)x rotação (rpm)")

    plt.figure()
    plt.plot(rotacao,Rendimento)
    plt.title("Rendimento x rotação (rpm)")

    plt.figure()
    plt.plot(rotacao,Tele)
    plt.title("Torque motor x rotação (rpm)")

    plt.figure()
    plt.plot(rotacao,Pmec)
    plt.title("Potencia Convertida x rotação (rpm)")

    plt.show()
if __name__=="__main__":
    main()
