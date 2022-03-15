import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares,leastsq

def FEMa(Iexc,c):
    VMax,VRem,K1,K2 = c
    return VMax*np.arctan(K1*Iexc)/(np.pi/2) + VRem*np.exp(-K2*Iexc)

def erro_FEMa(c,Iexc,Vexp):
    return FEMa(Iexc,c) - Vexp

def main():
    print('+======= Experimento 2: Curva de magnetização')
    df = pd.read_csv("data/flux.csv",skiprows=2)

    print("\n+== Dados:")
    print(df)

    # Dados Experimentais
    x = df["Ic_MCC"]
    y = df["Va_MCC"]
    e = df["E0_MCC"]

    # Definicao de constantes iniciais
    Vnom = 180
    IVn2 = df["Ic_MCC"][6]
    VMax0   = 1.1*Vnom
    VRem0   = y[0] #V(Iexc = 0)
    K10     = (4/np.pi) * 1/(IVn2)
    K20     = 3/.102
    
    c0 =  [VMax0,VRem0,K10,K20]

    mksize = 5 # Tamanho do ponto nos gráficos

    # Plot E0xI
    plt.figure()
    plt.plot(x, e     , marker = "o", markersize=mksize, label="$E_0$")
    plt.legend()
    plt.xlabel("Corrente [A]")
    plt.ylabel("Tensão E0 [V]")
    plt.grid()
    plt.savefig("img/E0.png")

    Va_0 = FEMa(x, c = c0)
    res,_ = leastsq(erro_FEMa, x0=c0, args=(x,y))
    Va_lsq = FEMa(x, res)

    ey0 = erro_FEMa(c0,Vexp=e, Iexc=x)
    
    print("\n+== Parâmetros\n  Vmax = {:.3f}\n  VRem = {:.3f}\n  K1   = {:.3f}\n  K2   = {:.3f}".format(*c0))
    print("\n+== Parâmetros(Least Squares Opt)\n  Vmax = {:.3f}\n  VRem = {:.3f}\n  K1   = {:.3f}\n  K2   = {:.3f}".format(*res))
    print("\n+== Pos escorv\n  Va_exp = {:.3f}\n  Va(Iexp={:.3f}) = {:.3f} ({:.3f} %)".format(df["Va_MCC"][6],.102,FEMa(.102,res),abs(1-(FEMa(.102,res)/df["Va_MCC"][6]))*100))

    ey0_res = erro_FEMa(res,Vexp=y, Iexc=x)

    # Plot Tensões
    plt.figure()
    plt.plot(x, y      , marker = "o", markersize=mksize, label="Experimental "+r"($V_a$)")
    plt.plot(x, Va_0   , marker = "o", markersize=mksize, label=r"$V_a(I_{exc})$")
    plt.plot(x, Va_lsq , marker = "o", markersize=mksize, label=r"$V_{aLSq}(I_{exc})$")
    plt.legend()
    plt.xlabel("Corrente [A]")
    plt.ylabel("Tensão [V]")
    plt.grid()
    plt.savefig("img/V_interp.png")

    # Plot Erro
    # plt.subplot(212)
    plt.figure()
    plt.title("Erro")
    plt.plot(x, ey0     , marker = "o", markersize=mksize, label="Erro "+r"$V_a$")
    plt.plot(x, ey0_res , marker = "o", markersize=mksize, label="Erro "+r"$V_{aLSq}$")
    plt.legend()
    plt.xlabel("Corrente [A]")
    plt.ylabel("Erro absoluto [V]")
    plt.grid()
    plt.savefig("img/V_interp_err.png")
    # plt.subplot(212)

    # Plot Erro percentual
    plt.figure()
    plt.title("Erro")
    plt.plot(x, (ey0/y)*100     , marker = "o", markersize=mksize, label="Erro "+r"$V_a$")
    plt.plot(x, (ey0_res/y)*100 , marker = "o", markersize=mksize, label="Erro "+r"$V_{aLSq}$")
    plt.legend()
    plt.xlabel("Corrente [A]")
    plt.ylabel("Erro percentual [%]")
    plt.grid()
    plt.savefig("img/V_interp_err_perc.png")
    # plt.show()


if __name__=="__main__":
    main()
