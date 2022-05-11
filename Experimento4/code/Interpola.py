import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares,leastsq
import os
import argparse
import tabulate

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

def experimento_4(Xm, Rf, R1, R2, X1, X2=None, base="experimento", save=False ,debug=False):
    # entre aqui os parâmetros da sua máquina
    Vfase = 220/np.sqrt(3) # esta é a tensão de fase do modelo.
    frequencia = 60#
    p = 2 # numero de par de polos
    We = 2*np.pi*frequencia # frequência angular das grandezas elétricas rad/s
    Ws = We/p # esta é a rotação síncrona da máquina em rad/s
    RPM=frequencia*60/p  # esta é a rotação síncrona da máquina, em RPM

    # Parâmetros a serem introduzidos
    if X2==None: X2=X1
    Par=[p, Xm, Rf, R1, R2, X1, X2]


    if debug: print("Gerando escorregamento")
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

    if debug: print("Gerando curva")
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
    if save:
        # # graficos em funcao da velocidade angular(pu)
        plt.figure()
        plt.title("Ia(A) x rotação")
        plt.subplot(2,1,1)
        plt.plot(rotacao,ModIa)
        plt.grid()
        plt.xlabel('Rotação [RPM]')
        plt.ylabel("Corrente de alimentação [A]")
        plt.subplot(2,1,2)
        plt.plot(rotacao,cosfi)
        plt.title("cos(fi)x rotação")
        plt.grid()
        plt.xlabel('Rotação [RPM]')
        plt.ylabel("cos(fi)")
        plt.savefig(f"img/{base}_Ia.png")

        plt.figure()
        plt.plot(rotacao,Rendimento)
        plt.title("Rendimento x rotação")
        plt.grid()
        plt.xlabel('Rotação [RPM]')
        plt.ylabel("Rendimento")
        plt.savefig(f"img/{base}_Rendimento.png")

        plt.figure()
        plt.plot(rotacao,Tele)
        plt.title("Torque motor x rotação")
        plt.grid()
        plt.xlabel('Rotação [RPM]')
        plt.ylabel("Torque [Nm]")
        plt.savefig(f"img/{base}_Torque.png")

        plt.figure()
        plt.plot(rotacao,Pmec)
        plt.title("Potencia Convertida x rotação")
        plt.grid()
        plt.xlabel('Rotação [RPM]')
        plt.ylabel("Potencia Convertida [W]")
        plt.savefig(f"img/{base}_Potencia.png")

    return rotacao, ModIa, cosfi, Rendimento, Tele, Pmec

def n_perc(Pconv,V,I,fp):
    return 100*Pconv/(np.sqrt(3)*V*I*fp)
    
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
    print("+ Calculando para método Aurio")
    Xm=10.03
    Rf=58.04
    R1=2.7
    R2=3.5
    X1=1.95
    rotacao_aurio, ModIa_aurio, cosfi_aurio, Rendimento_aurio, Tele_aurio, Pmec_aurio = experimento_4(Xm, Rf, R1, R2, X1, base="aurio", save=True) 

    print("+ Calculando para método IEEE")
    Xm=15.98
    Rf=58.04
    R1=2.7
    R2=3.5
    X1=1.95
    rotacao_ieee, ModIa_ieee, cosfi_ieee, Rendimento_ieee, Tele_ieee, Pmec_ieee = experimento_4(Xm, Rf, R1, R2, X1, base="ieee", save=True)    

    base = "aurio_ieee"
    # # graficos em funcao da velocidade angular(pu)
    plt.figure(figsize=[6.4,9])
    plt.subplot(2,1,1)
    plt.title("Ia(A) x rotação")
    plt.plot(rotacao_aurio,ModIa_aurio, label="Aurio")
    plt.plot(rotacao_ieee,ModIa_ieee, label="IEEE")
    plt.grid()
    plt.xlabel('Rotação [RPM]')
    plt.ylabel("Corrente de alimentação [A]")
    plt.subplot(2,1,2)
    plt.plot(rotacao_aurio,cosfi_aurio, label="Aurio")
    plt.plot(rotacao_ieee,cosfi_ieee, label="IEEE")
    plt.title("cos(fi)x rotação")
    plt.grid()
    plt.xlabel('Rotação [RPM]')
    plt.ylabel("cos(fi)")
    plt.legend()
    plt.savefig(f"img/{base}_Ia.png")

    plt.figure()
    plt.plot(rotacao_aurio,Rendimento_aurio, label="Aurio")
    plt.plot(rotacao_ieee,Rendimento_ieee, label="IEEE")
    plt.title("Rendimento x rotação")
    plt.grid()
    plt.xlabel('Rotação [RPM]')
    plt.ylabel("Rendimento")
    plt.legend()
    plt.savefig(f"img/{base}_Rendimento.png")

    plt.figure()
    plt.plot(rotacao_aurio,Tele_aurio, label="Aurio")
    plt.plot(rotacao_ieee,Tele_ieee, label="IEEE")
    plt.title("Torque motor x rotação")
    plt.grid()
    plt.xlabel('Rotação [RPM]')
    plt.ylabel("Torque [Nm]")
    plt.legend()
    plt.savefig(f"img/{base}_Torque.png")

    plt.figure()
    plt.plot(rotacao_aurio,Pmec_aurio, label="Aurio")
    plt.plot(rotacao_ieee,Pmec_ieee, label="IEEE")
    plt.title("Potencia Convertida x rotação")
    plt.grid()
    plt.xlabel('Rotação [RPM]')
    plt.ylabel("Potencia Convertida [W]")
    plt.legend()
    plt.savefig(f"img/{base}_Potencia.png")


    print("+ Definindo pontos para 750W")
    eps = 2
    print(f"  Threshold de {eps} W.")

    result_aurio = []
    for rot,pot in zip(rotacao_aurio,Pmec_aurio):
        if abs(pot-750)<eps:
            result_aurio.append((rot,pot,pot-750))
    result_ieee = []
    for rot,pot in zip(rotacao_ieee,Pmec_ieee):
        if abs(pot-750)<eps:
            result_ieee.append((rot,pot,pot-750))
    
    print(f"  AURIO - Encontrados {len(result_aurio)} resultados:")
    for res in result_aurio:
        print(f"    {res}")

    print(f"  IEEE - Encontrados {len(result_ieee)} resultados:")
    for res in result_ieee:
        print(f"    {res}")

    rot_750_aurio, Pmec_750_aurio, _ = result_aurio[-1]
    rot_750_ieee, Pmec_750_ieee, _ = result_ieee[-1]
    
    torque_750_aurio = Tele_aurio[rotacao_aurio.index(rot_750_aurio)]
    torque_750_ieee = Tele_ieee[rotacao_ieee.index(rot_750_ieee)]
    
    Ia_750_aurio = ModIa_aurio[rotacao_aurio.index(rot_750_aurio)]
    Ia_750_ieee = ModIa_ieee[rotacao_ieee.index(rot_750_ieee)]

    fp_750_aurio = cosfi_aurio[rotacao_aurio.index(rot_750_aurio)]
    fp_750_ieee = cosfi_ieee[rotacao_ieee.index(rot_750_ieee)]
    
    n_perc_750_aurio = n_perc(Pconv = Pmec_750_aurio, V =  220/np.sqrt(3), I = Ia_750_aurio, fp = fp_750_aurio)
    n_perc_750_ieee  = n_perc(Pconv = Pmec_750_ieee,  V =  220/np.sqrt(3), I = Ia_750_ieee,  fp = fp_750_ieee)

    dict_750={"Aurio":{
                    "Torque": torque_750_aurio,
                    "Rotação": rot_750_aurio,
                    "FP": fp_750_aurio,
                    "n%" : n_perc_750_aurio

                },
                "IEEE":{
                    "Torque": torque_750_ieee,
                    "Rotação": rot_750_ieee,
                    "FP": fp_750_ieee,
                    "n%" : n_perc_750_ieee
                },}

    df = pd.DataFrame(dict_750).T
    print(tabulate.tabulate(df,headers=df.columns, tablefmt="psql"))
    

if __name__=="__main__":
    main()
