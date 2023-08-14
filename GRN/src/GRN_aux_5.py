import matplotlib.pyplot as plt
import copy
import random as r
import numpy as np
from scipy.integrate import odeint
import math
import matplotlib

arquivo = open("GRN5.txt", 'r')

x = []
A = []
B = []
C = []
D = []
E = []

for linha in arquivo:
    elementos = linha.split()
    x.append(float(elementos[0].strip()))
    A.append(float(elementos[1].strip()))
    B.append(float(elementos[2].strip()))
    C.append(float(elementos[3].strip()))
    D.append(float(elementos[4].strip()))
    E.append(float(elementos[5].strip()))

maximo_A = max(A)
maximo_B = max(B)
maximo_C = max(C)
maximo_D = max(D)
maximo_E = max(E)

maximos = []
maximos.append(maximo_A)
maximos.append(maximo_B)
maximos.append(maximo_C)
maximos.append(maximo_D)
maximos.append(maximo_E)

Y0 = []
Y0.append(A[0])
Y0.append(B[0])
Y0.append(C[0])
Y0.append(D[0])
Y0.append(E[0])

IND_SIZE = 15  # Tamanho do indivíduo (quantidade de coeficientes)
MIN_K = 0.1  # Menor valor que K pode assumir
MAX_K = 1  # Maior valor que K pode assumir
MIN_N = 1  # Menor valor que N pode assumir
MAX_N = 25  # Maior valor que N pode assumir
MIN_TAU = 0.1  # Menor valor que TAU pode assumir
MAX_TAU = 5  # Maior valor que TAU pode assumir
MIN_STRATEGY = 0.1  # Menor valor que a estratégia pode assumir
MAX_STRATEGY = 10  # Maior valor que a estratégia pode assumir
TAU_SIZE = 5
N_SIZE = 5
K_SIZE = 5
# INDIVIDUO: [tauA, tauB, tauC, tauD, tauE, tauF, tauG, tauH, tauI, tauJ, kAJ, kBE, kCB, kCF, kCA, kDF, kEJ, kFA, kGB, kGF, kGA, kHF, kIG, kIH, kJI, nAJ, nBE, nCB, nCF, nCA, nDF, nEJ, nFA, nGB, nGF, nGA, nHF, nIG, nIH, nJI]


# dobra_pontos = [0.0, 0.7347, 1.4694, 2.2041, 2.9388, 3.6735, 4.4082, 5.1429, 5.8776, 6.6123, 7.3469, 8.0816, 8.8163, 9.551, 10.2857, 11.0204, 11.7551, 12.4898, 13.2245, 13.9592, 14.6939, 15.4286, 16.1633, 16.898, 17.6327, 18.3674, 19.102, 19.8367, 20.5714, 21.3061, 22.0408, 22.7755, 23.5102, 24.2449, 24.9796, 25.7143, 26.449, 27.1837, 27.9184, 28.6531, 29.3878, 30.1225, 30.8571, 31.5918, 32.3265, 33.0612, 33.7959, 34.5306, 35.2653, 36.0, 36.7347, 37.4694, 38.2041, 38.9388, 39.6735, 40.4082, 41.1429, 41.8776, 42.6122, 43.3469, 44.0816, 44.8163, 45.551, 46.2857, 47.0204, 47.7551, 48.4898, 49.2245, 49.9592, 50.6939, 51.4286, 52.1633, 52.898, 53.6327, 54.3673, 55.102, 55.8367, 56.5714, 57.3061, 58.0408, 58.7755, 59.5102, 60.2449, 60.9796, 61.7143, 62.449, 63.1837, 63.9184, 64.6531, 65.3878, 66.1224, 66.8571, 67.5918, 68.3265, 69.0612, 69.7959, 70.5306, 71.2653, 72.0, 72.7347]
dobra_pontos = copy.deepcopy(x)


def twoBody(y, t, tauA, kA, nA, tauB, kB, nB, tauC, kC, nC, tauD, kD, nD, tauE, kE, nE):
    ydot = np.empty((5,))

    ydot[0] = ((1 - (pow((y[4] / maximo_E), nA)) / (pow((y[4] / maximo_E), nA) + pow(kA, nA))) - (
                y[0] / maximo_A)) / tauA

    ydot[1] = (((pow((y[0] / maximo_A), nB)) / (pow((y[0] / maximo_A), nB) + pow(kB, nB))) - (y[1] / maximo_B)) / tauB

    ydot[2] = (((pow((y[1] / maximo_B), nC)) / (pow((y[1] / maximo_B), nC) + pow(kC, nC))) - (y[2] / maximo_C)) / tauC

    ydot[3] = (((pow((y[2] / maximo_C), nD)) / (pow((y[2] / maximo_C), nD) + pow(kD, nD))) - (y[3] / maximo_D)) / tauD

    ydot[4] = (((pow((y[3] / maximo_D), nE)) / (pow((y[3] / maximo_D), nE) + pow(kE, nE))) - (y[4] / maximo_E)) / tauE

    return ydot


def twoBodyNew(y, t, tauA, kA, nA, tauB, kB, nB, tauC, kC, nC, tauD, kD, nD, tauE, kEB, kED, kEE, nEB, nED, nEE):
    ydot = np.empty((5,))

    ydot[0] = ((1 - (pow((y[4] / maximo_E), nA) / (pow((y[4] / maximo_E), nA) + pow(kA, nA)))) - (
                y[0] / maximo_A)) / tauA

    ydot[1] = (((pow((y[0] / maximo_A), nB)) / (pow((y[0] / maximo_A), nB) + pow(kB, nB))) - (y[1] / maximo_B)) / tauB

    ydot[2] = (((pow((y[1] / maximo_B), nC)) / (pow((y[1] / maximo_B), nC) + pow(kC, nC))) - (y[2] / maximo_C)) / tauC

    ydot[3] = (((pow((y[2] / maximo_C), nD)) / (pow((y[2] / maximo_C), nD) + pow(kD, nD))) - (y[3] / maximo_D)) / tauD

    ydot[4] = ((((pow(y[1] / maximo_B, nEB) / (pow(y[1] / maximo_B, nEB) + pow(kEB, nEB))) * (
                pow(y[3] / maximo_D, nED) / (pow(y[3] / maximo_D, nED) + pow(kED, nED)))) + (
                            (pow(y[3] / maximo_D, nED) / (pow(y[3] / maximo_D, nED) + pow(kED, nED))) * (
                                pow(y[4] / maximo_E, nEE) / (pow(y[4] / maximo_E, nEE) + pow(kEE, nEE))))) - (
                           y[4] / maximo_E)) / tauE

    return ydot

def organiza_pontos(solucao):
    pA = []
    pB = []
    pC = []
    pD = []
    pE = []
    for pontos in range(len(solucao)):
        if pontos % 2 == 0 or pontos % 2 == 1:
            pA.append(solucao[pontos][0])
            pB.append(solucao[pontos][1])
            pC.append(solucao[pontos][2])
            pD.append(solucao[pontos][3])
            pE.append(solucao[pontos][4])
    return pA, pB, pC, pD, pE


def calcula_diferenca(pA, pB, pC, pD, pE):  # , pF, pG, pH, pI, pJ):
    difA = 0
    difB = 0
    difC = 0
    difD = 0
    difE = 0


    for elemento in range(len(pA)):
        dif = abs(A[elemento] - pA[elemento])
        difA += dif
    for elemento in range(len(pB)):
        dif = abs(B[elemento] - pB[elemento])
        difB += dif
    for elemento in range(len(pC)):
        dif = abs(C[elemento] - pC[elemento])
        difC += dif
    for elemento in range(len(pD)):
        dif = abs(D[elemento] - pD[elemento])
        difD += dif
    for elemento in range(len(pE)):
        dif = abs(E[elemento] - pE[elemento])
        difE += dif

    difTotal = difA + difB + difC + difD + difE  # + difF + difG + difH + difI + difJ
    # difTotal = abs(difTotal)
    return difTotal



def evaluation_5(ind_atual):
    #ind_atual = POPULACAO[qtd_progenitores]
    tauA = ind_atual[0]
    tauB = ind_atual[1]
    tauC = ind_atual[2]
    tauD = ind_atual[3]
    tauE = ind_atual[4]
    kA = ind_atual[5]
    kB = ind_atual[6]
    kC = ind_atual[7]
    kD = ind_atual[8]
    kE = ind_atual[9]
    nA = ind_atual[10]
    nB = ind_atual[11]
    nC = ind_atual[12]
    nD = ind_atual[13]
    nE = ind_atual[14]


    solution = odeint(twoBody, Y0, dobra_pontos, args=(
        tauA, kA, int(nA), tauB, kB, int(nB), tauC, kC, int(nC), tauD, kD, int(nD), tauE, kE, int(nE),))
    pA, pB, pC, pD, pE = organiza_pontos(solution)
    return calcula_diferenca(pA, pB, pC, pD, pE)


def get_result_static():
     #ind_atual = POPULACAO[qtd_progenitores]
    tauA = 1.25
    tauB = 4
    tauC = 1.02
    tauD = 1.57
    tauE = 3.43
    kA = 0.72
    kB = 0.5
    kC = 0.45
    kD = 0.51
    kE = 0.52
    nA = 13
    nB = 4
    nC = 3
    nD = 4
    nE = 16


    solution = odeint(twoBody, Y0, dobra_pontos, args=(
        tauA, kA, int(nA), tauB, kB, int(nB), tauC, kC, int(nC), tauD, kD, int(nD), tauE, kE, int(nE),))
    pA, pB, pC, pD, pE = organiza_pontos(solution)
    #return calcula_diferenca(pA, pB, pC, pD, pE)

    return pA, pB, pC, pD, pE

def get_results(ind_atual):
    #ind_atual = POPULACAO[qtd_progenitores]
    tauA = ind_atual[0]
    tauB = ind_atual[1]
    tauC = ind_atual[2]
    tauD = ind_atual[3]
    tauE = ind_atual[4]
    kA = ind_atual[5]
    kB = ind_atual[6]
    kC = ind_atual[7]
    kD = ind_atual[8]
    kE = ind_atual[9]
    nA = ind_atual[10]
    nB = ind_atual[11]
    nC = ind_atual[12]
    nD = ind_atual[13]
    nE = ind_atual[14]


    solution = odeint(twoBody, Y0, dobra_pontos, args=(
        tauA, kA, int(nA), tauB, kB, int(nB), tauC, kC, int(nC), tauD, kD, int(nD), tauE, kE, int(nE),))
    pA, pB, pC, pD, pE = organiza_pontos(solution)

    print(calcula_diferenca(pA, pB, pC, pD, pE))
    return pA, pB, pC, pD, pE

def get_results_new(ind_atual):
   # ind_atual = [1.2163355099083872, 1.1264485098219865, 2.973714367061704, 2.952143123315177, 2.998260518457365, 0.5687249950503857, 0.4580723119903261, 0.46214892372246563, 0.6182568295500336, 0.5213082492659304, 0.7708877748759901, 0.1497642024548283, 4.254757908429968, 3.759370669969996, 4.784173526119725, 10.935884810737809, 24.595975874929724, 2.8109199678182635, 4.922623602327875, 1.804297289687443, 0.6961641316460799, 1.1805067448542073, 3.769380877770944, 0.6268588518301711, 0.7945926074279098, 0.8665901143646684, 1.1631405647512596, 2.3141178250393146, 1.2841062086785697, 0.7091739090057955, 2.245324305294896, 0.8208096283146853, 0.6911548119817139, 0.7874359961268611, 1.101498884676361, 0.7742237807425528, 0.9816157239798934, 70.97756457350062]

    tauA = ind_atual[0]
    tauB = ind_atual[1]
    tauC = ind_atual[2]
    tauD = ind_atual[3]
    tauE = ind_atual[4]
    kA = ind_atual[5]
    kB = ind_atual[6]
    kC = ind_atual[7]
    kD = ind_atual[8]
    kEB = ind_atual[9]
    kED = ind_atual[10]
    kEE = ind_atual[11]
    nA = ind_atual[12]
    nB = ind_atual[13]
    nC = ind_atual[14]
    nD = ind_atual[15]
    nEB = ind_atual[16]
    nED = ind_atual[17]
    nEE = ind_atual[18]


    solution = odeint(twoBodyNew, Y0, dobra_pontos, args=(tauA, kA, int(nA), tauB, kB, int(nB), tauC, kC, int(nC), tauD, kD, int(nD), tauE, kEB, kED, kEE, int(nEB), int(nED), int(nEE),))
    pA, pB, pC, pD, pE = organiza_pontos(solution)

    print(calcula_diferenca(pA, pB, pC, pD, pE))
    return pA, pB, pC, pD, pE

def diferenca_teste():
    tauA = 1.25
    tauB = 4
    tauC = 1.02
    tauD = 1.57
    tauE = 3.43
    kA = 0.72
    kB = 0.5
    kC = 0.45
    kD = 0.51
    kE = 0.52
    nA = 13
    nB = 4
    nC = 3
    nD = 4
    nE = 16

    solution = odeint(twoBody, Y0, dobra_pontos, args=(
        tauA, kA, int(nA), tauB, kB, int(nB), tauC, kC, int(nC), tauD, kD, int(nD), tauE, kE, int(nE),))
    pA, pB, pC, pD, pE = organiza_pontos(solution)

    print(solution)

    return calcula_diferenca(pA, pB, pC, pD, pE)