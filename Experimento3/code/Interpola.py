import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares,leastsq
import os
import argparse

def Eab(Iexc,constants):
    Vab_max,Vab_rem,k1,k2 = constants
    return Vab_max*np.arctan(k1*Iexc)/(np.pi/2) + Vab_rem*np.exp(-k2*Iexc)

def erro_Eab(constants,Iexc,Vexp):
    return Eab(Iexc,constants) - Vexp

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
    print('+======= Experimento 3: ')
    df = pd.read_csv(data_file, decimal=",")

    print("\n+== Dados:")
    print(df)

    # Dados Experimentais
    x = df["Icampo (A)"]#   
    y = df["Va de linha medido(V)"]
    e = df["Va de fase calculado (V)"]
    
    # # Definicao de constantes iniciais
    Vnom = 180
    IVn2 = x[6]
    VMax0   = Vnom
    VRem0   = y[0] #V(Iexc = 0)
    K10     = 1
    K20     = 5
    
    c0 =  [VMax0,VRem0,K10,K20]

    Va_0 = Eab(x, constants = c0)
    res,_ = leastsq(erro_Eab, x0=c0, args=(x,y))
    Va_lsq = Eab(x, res)

    ey0 = erro_Eab(c0,Vexp=e, Iexc=x)
    
    print("\n+== Parâmetros\n  Vmax = {:.3f}\n  VRem = {:.3f}\n  K1   = {:.3f}\n  K2   = {:.3f}".format(*c0))
    print("\n+== Parâmetros(Least Squares Opt)\n  Vmax = {:.3f}\n  VRem = {:.3f}\n  K1   = {:.3f}\n  K2   = {:.3f}".format(*res))

    ey0_res = erro_Eab(res,Vexp=y, Iexc=x)

    mksize = 5 # Tamanho do ponto nos gráficos
    plt.figure()
    plt.plot(x, y, marker = "o", linestyle="", markersize=mksize, label="Medido")
    plt.plot(x, e, marker = "o", linestyle="", markersize=mksize, label="Calculado")
    plt.plot(x, Va_0, marker = "o", linestyle="", markersize=mksize, label="Calculado Py")
    plt.legend()
    plt.xlabel("Corrente de campo [A]")
    plt.ylabel("Tensão de armadura [V]")
    plt.grid()
    if(args.save): plt.savefig(os.path.join(img_path,"E0.png"))
    
    # Plot Tensões
    plt.figure()
    plt.plot(x, y      , marker = "o", markersize=mksize, label="Experimental "+r"($V_a$)")
    plt.plot(x, Va_0   , marker = "o", markersize=mksize, label=r"$V_a(I_{exc})$")
    plt.plot(x, Va_lsq , marker = "o", markersize=mksize, label=r"$V_{aLSq}(I_{exc})$")
    plt.legend()
    plt.xlabel("Corrente [A]")
    plt.ylabel("Tensão [V]")
    plt.grid()
    if(args.save): plt.savefig(os.path.join(img_path,"V_interp.png"))

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
    if(args.save): plt.savefig(os.path.join(img_path,"V_interp_err.png"))
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
    if(args.save): plt.savefig(os.path.join(img_path,"V_interp_err_perc.png"))
    """
    """
    if not args.save: plt.show()

if __name__=="__main__":
    main()
