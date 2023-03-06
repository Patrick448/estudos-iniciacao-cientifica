import matplotlib.pyplot as plt
import copy
import random as r
import numpy as np
from scipy.integrate import odeint
import math
import matplotlib


#Leitura do arquivo com os pontos de expressão gênica
arquivo = open("GRN10.txt", 'r')

x = []
A = []
B = []
C = []
D = []
E = []
F = []
G = []
H = []
I = []
J = []

for linha in arquivo:
    elementos = linha.split()
    x.append(float(elementos[0].strip()))
    A.append(float(elementos[1].strip()))
    B.append(float(elementos[2].strip()))
    C.append(float(elementos[3].strip()))
    D.append(float(elementos[4].strip()))
    E.append(float(elementos[5].strip()))
    F.append(float(elementos[6].strip()))
    G.append(float(elementos[7].strip()))
    H.append(float(elementos[8].strip()))
    I.append(float(elementos[9].strip()))
    J.append(float(elementos[10].strip()))


#Valores máximos da expressão gênica de cada espécie
maximo_A = max(A)
maximo_B = max(B)
maximo_C = max(C)
maximo_D = max(D)
maximo_E = max(E)
maximo_F = max(F)
maximo_G = max(G)
maximo_H = max(H)
maximo_I = max(I)
maximo_J = max(J)


#Condições iniciais da EDO
Y0 = []
Y0.append(A[0])
Y0.append(B[0])
Y0.append(C[0])
Y0.append(D[0])
Y0.append(E[0])
Y0.append(F[0])
Y0.append(G[0])
Y0.append(H[0])
Y0.append(I[0])
Y0.append(J[0])

#Configuração do Indivíduo e range de parâmetros
IND_SIZE = 40#Tamanho do indivíduo (quantidade de coeficientes)
MIN_K = 0.1#Menor valor que K pode assumir
MAX_K = 1#Maior valor que K pode assumir
MIN_N = 1#Menor valor que N pode assumir
MAX_N = 25#Maior valor que N pode assumir
MIN_TAU = 0.1#Menor valor que TAU pode assumir
MAX_TAU = 5#Maior valor que TAU pode assumir
MIN_STRATEGY = 0.1 #Menor valor que a estratégia pode assumir
MAX_STRATEGY = 10#Maior valor que a estratégia pode assumir
TAU_SIZE = 10
N_SIZE = 15
K_SIZE = 15
#COMPOSIÇÃO DO INDIVIDUO: [tauA, tauB, tauC, tauD, tauE, tauF, tauG, tauH, tauI, tauJ, kAJ, kBE, kCB, kCF, kCA, kDF, kEJ, kFA, kGB, kGF, kGA, kHF, kIG, kIH, kJI, nAJ, nBE, nCB, nCF, nCA, nDF, nEJ, nFA, nGB, nGF, nGA, nHF, nIG, nIH, nJI]


dobra_pontos = copy.deepcopy(x)
POPULACAO = []
APTIDAO = []


LAMBDA_FILHOS = []
APTIDAO_FILHOS = []


#Função de integração numérica - cada ydot é uma EDO de uma espécie 
def twoBody(y, t, tauA, kAJ, nAJ, tauB, kBE, nBE, tauC, kCB, nCB, kCF, nCF, kCA, nCA, tauD, kDF, nDF, tauE, kEJ, nEJ, tauF, kFA, nFA, tauG, kGB, nGB, kGF, nGF, kGA, nGA, tauH, kHF, nHF, tauI, kIG, nIG, kIH, nIH, tauJ, kJI, nJI):
    ydot = np.empty((10,))

    ydot[0] = ((1-pow(y[9]/maximo_J,nAJ)/(pow(y[9]/maximo_J,nAJ)+pow(kAJ,nAJ)))-(y[0]/maximo_A)) / tauA

    ydot[1] = (pow(y[4]/maximo_E,nBE)/(pow(y[4]/maximo_E,nBE)+pow(kBE,nBE))-(y[1]/maximo_B)) / tauB

    ydot[2] = (pow(y[1]/maximo_B,nCB)/(pow(y[1]/maximo_B,nCB)+pow(kCB,nCB))*(1-pow(y[5]/maximo_F,nCF)/(pow(y[5]/maximo_F,nCF)+pow(kCF,nCF)))*(1-pow(y[0]/maximo_A,nCA)/(pow(y[0]/maximo_A,nCA)+pow(kCA,nCA)))+(1-pow(y[1]/maximo_B,nCB)/(pow(y[1]/maximo_B,nCB)+pow(kCB,nCB)))*pow(y[5]/maximo_F,nCF)/(pow(y[5]/maximo_F,nCF)+pow(kCF,nCF))*(1-pow(y[0]/maximo_A,nCA)/(pow(y[0]/maximo_A,nCA)+pow(kCA,nCA)))+(1-pow(y[1]/maximo_B,nCB)/(pow(y[1]/maximo_B,nCB)+pow(kCB,nCB)))*(1-pow(y[5]/maximo_F,nCF)/(pow(y[5]/maximo_F,nCF)+pow(kCF,nCF)))*pow(y[0]/maximo_A,nCA)/(pow(y[0]/maximo_A,nCA)+pow(kCA,nCA))+pow(y[1]/maximo_B,nCB)/(pow(y[1]/maximo_B,nCB)+pow(kCB,nCB))*(1-pow(y[5]/maximo_F,nCF)/(pow(y[5]/maximo_F,nCF)+pow(kCF,nCF)))*pow(y[0]/maximo_A,nCA)/(pow(y[0]/maximo_A,nCA)+pow(kCA,nCA))+(1-pow(y[1]/maximo_B,nCB)/(pow(y[1]/maximo_B,nCB)+pow(kCB,nCB)))*pow(y[5]/maximo_F,nCF)/(pow(y[5]/maximo_F,nCF)+pow(kCF,nCF))*pow(y[0]/maximo_A,nCA)/(pow(y[0]/maximo_A,nCA)+pow(kCA,nCA))+pow(y[1]/maximo_B,nCB)/(pow(y[1]/maximo_B,nCB)+pow(kCB,nCB))*pow(y[5]/maximo_F,nCF)/(pow(y[5]/maximo_F,nCF)+pow(kCF,nCF))*pow(y[0]/maximo_A,nCA)/(pow(y[0]/maximo_A,nCA)+pow(kCA,nCA))-(y[2]/maximo_C)) / tauC

    ydot[3] = (pow(y[5]/maximo_F,nDF)/(pow(y[4]/maximo_E,nDF)+pow(kDF,nDF))-(y[3]/maximo_D)) / tauD

    ydot[4] = (1-pow(y[9]/maximo_J,nEJ)/(pow(y[9]/maximo_J,nEJ)+pow(kEJ,nEJ))-(y[4]/maximo_E)) / tauE

    ydot[5] = (pow(y[0]/maximo_A,nFA)/(pow(y[0]/maximo_A,nFA)+pow(kFA,nFA))-(y[5]/maximo_F)) / tauF

    ydot[6] = (pow(y[1]/maximo_B,nGB)/(pow(y[1]/maximo_B,nGB)+pow(kGB,nGB))*(1-pow(y[5]/maximo_F,nGF)/(pow(y[5]/maximo_F,nGF)+pow(kGF,nGF)))*(1-pow(y[0]/maximo_A,nGA)/(pow(y[0]/maximo_A,nGA)+pow(kGA,nGA)))+(1-pow(y[1]/maximo_B,nGB)/(pow(y[1]/maximo_B,nGB)+pow(kGB,nGB)))*pow(y[5]/maximo_F,nGF)/(pow(y[5]/maximo_F,nGF)+pow(kGF,nGF))*(1-pow(y[0]/maximo_A,nGA)/(pow(y[0]/maximo_A,nGA)+pow(kGA,nGA)))+(1-pow(y[1]/maximo_B,nGB)/(pow(y[1]/maximo_B,nGB)+pow(kGB,nGB)))*(1-pow(y[5]/maximo_F,nGF)/(pow(y[5]/maximo_F,nGF)+pow(kGF,nGF)))*pow(y[0]/maximo_A,nGA)/(pow(y[0]/maximo_A,nGA)+pow(kGA,nGA))+pow(y[1]/maximo_B,nGB)/(pow(y[1]/maximo_B,nGB)+pow(kGB,nGB))*(1-pow(y[5]/maximo_F,nGF)/(pow(y[5]/maximo_F,nGF)+pow(kGF,nGF)))*pow(y[0]/maximo_A,nGA)/(pow(y[0]/maximo_A,nGA)+pow(kGA,nGA))+(1-pow(y[1]/maximo_B,nGB)/(pow(y[1]/maximo_B,nGB)+pow(kGB,nGB)))*pow(y[5]/maximo_F,nGF)/(pow(y[5]/maximo_F,nGF)+pow(kGF,nGF))*pow(y[0]/maximo_A,nGA)/(pow(y[0]/maximo_A,nGA)+pow(kGA,nGA))+pow(y[1]/maximo_B,nGB)/(pow(y[1]/maximo_B,nGB)+pow(kGB,nGB))*pow(y[5]/maximo_F,nGF)/(pow(y[5]/maximo_F,nGF)+pow(kGF,nGF))*pow(y[0]/maximo_A,nGA)/(pow(y[0]/maximo_A,nGA)+pow(kGA,nGA))-(y[6]/maximo_G)) / tauG

    ydot[7] = (pow(y[5]/maximo_F,nHF)/(pow(y[5]/maximo_F,nHF)+pow(kHF,nHF))-(y[7]/maximo_H)) / tauH

    ydot[8] = (pow(y[6]/maximo_G,nIG)/(pow(y[6]/maximo_G,nIG)+pow(kIG,nIG))*pow(y[7]/maximo_H,nIH)/(pow(y[7]/maximo_H,nIH)+pow(kIH,nIH))-(y[8]/maximo_I)) / tauI

    ydot[9] = (pow(y[8]/maximo_I,nJI)/(pow(y[8]/maximo_I,nJI)+pow(kJI,nJI))-(y[9]/maximo_J)) / tauJ

    return ydot

def organiza_pontos(solucao):
    pA = []
    pB = []
    pC = []
    pD = []
    pE = []
    pF = []
    pG = []
    pH = []
    pI = []
    pJ = []
    for pontos in range(len(solucao)):
        if pontos % 2 == 0 or pontos % 2 == 1:
            pA.append(solucao[pontos][0])
            pB.append(solucao[pontos][1])
            pC.append(solucao[pontos][2])
            pD.append(solucao[pontos][3])
            pE.append(solucao[pontos][4])
            pF.append(solucao[pontos][5])
            pG.append(solucao[pontos][6])
            pH.append(solucao[pontos][7])
            pI.append(solucao[pontos][8])
            pJ.append(solucao[pontos][9])
    return pA, pB, pC, pD, pE, pF, pG, pH, pI, pJ



def calcula_diferenca(pA, pB, pC, pD, pE, pF, pG, pH, pI, pJ):
    difA = 0
    difB = 0
    difC = 0
    difD = 0
    difE = 0
    difF = 0
    difG = 0
    difH = 0
    difI = 0
    difJ = 0
    pAl = []
    pBl = []
    pCl = []
    pDl = []
    pEl = []
    pFl = []
    pGl = []
    pHl = []
    pIl = []
    pJl = []
     
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
    for elemento in range(len(pF)):
        dif = abs(F[elemento] - pF[elemento])
        difF += dif
    for elemento in range(len(pG)):
        dif = abs(G[elemento] - pG[elemento])
        difG += dif
    for elemento in range(len(pH)):
        dif = abs(H[elemento] - pH[elemento])
        difH += dif
    for elemento in range(len(pI)):
        dif = abs(I[elemento] - pI[elemento])
        difI += dif
    for elemento in range(len(pJ)):
        dif = abs(J[elemento] - pJ[elemento])
        difJ += dif        
    
    difTotal = difA + difB + difC + difD + difE + difF + difG + difH + difI + difJ

    return difTotal




def evaluation(ind_atual):
    tauA = ind_atual[0]
    tauB = ind_atual[1]
    tauC = ind_atual[2]
    tauD = ind_atual[3]
    tauE = ind_atual[4]
    tauF = ind_atual[5]
    tauG = ind_atual[6]
    tauH = ind_atual[7]
    tauI = ind_atual[8]
    tauJ = ind_atual[9]
    kAJ = ind_atual[10]
    kBE = ind_atual[11]
    kCB = ind_atual[12]
    kCF = ind_atual[13]
    kCA = ind_atual[14]
    kDF = ind_atual[15]
    kEJ = ind_atual[16]
    kFA = ind_atual[17]
    kGB = ind_atual[18]
    kGF = ind_atual[19]
    kGA = ind_atual[20]
    kHF = ind_atual[21]
    kIG = ind_atual[22]
    kIH = ind_atual[23]
    kJI = ind_atual[24]
    nAJ = ind_atual[25]
    nBE = ind_atual[26]
    nCB = ind_atual[27]
    nCF = ind_atual[28]
    nCA = ind_atual[29]
    nDF = ind_atual[30]
    nEJ = ind_atual[31]
    nFA = ind_atual[32]
    nGB = ind_atual[33]
    nGF = ind_atual[34]
    nGA = ind_atual[35]
    nHF = ind_atual[36]
    nIG = ind_atual[37]
    nIH = ind_atual[38]
    nJI = ind_atual[39]

    solution = odeint(twoBody, Y0, dobra_pontos, args=(
    tauA, kAJ, int(nAJ), tauB, kBE, int(nBE), tauC, kCB, int(nCB), kCF, int(nCF), kCA, int(nCA), tauD, kDF, int(nDF),
    tauE, kEJ, int(nEJ), tauF, kFA, int(nFA), tauG, kGB, int(nGB), kGF, int(nGF), kGA, int(nGA), tauH, kHF, int(nHF),
    tauI, kIG, int(nIG), kIH, int(nIH), tauJ, kJI, int(nJI),))
    pA, pB, pC, pD, pE, pF, pG, pH, pI, pJ = organiza_pontos(solution)
    return calcula_diferenca(pA, pB, pC, pD, pE, pF, pG, pH, pI, pJ)

def diferenca_teste():

    ind_atual = [1.73,2,0.81,0.11, 1.23, 1.78, 1.14, 1.04, 3.47, 3.21,
                    0.45, 0.56, 0.99, 0.77, 0.71, 0.66, 0.46, 0.48, 0.66, 0.99, 0.85, 0.61, 0.55, 0.46, 0.17,
                    20, 9, 24, 12, 2, 2, 6, 4, 7, 24, 2, 7, 21, 20, 3 ]
    tauA = ind_atual[0]
    tauB = ind_atual[1]
    tauC = ind_atual[2]
    tauD = ind_atual[3]
    tauE = ind_atual[4]
    tauF = ind_atual[5]
    tauG = ind_atual[6]
    tauH = ind_atual[7]
    tauI = ind_atual[8]
    tauJ = ind_atual[9]
    kAJ = ind_atual[10]
    kBE = ind_atual[11]
    kCB = ind_atual[12]
    kCF = ind_atual[13]
    kCA = ind_atual[14]
    kDF = ind_atual[15]
    kEJ = ind_atual[16]
    kFA = ind_atual[17]
    kGB = ind_atual[18]
    kGF = ind_atual[19]
    kGA = ind_atual[20]
    kHF = ind_atual[21]
    kIG = ind_atual[22]
    kIH = ind_atual[23]
    kJI = ind_atual[24]
    nAJ = ind_atual[25]
    nBE = ind_atual[26]
    nCB = ind_atual[27]
    nCF = ind_atual[28]
    nCA = ind_atual[29]
    nDF = ind_atual[30]
    nEJ = ind_atual[31]
    nFA = ind_atual[32]
    nGB = ind_atual[33]
    nGF = ind_atual[34]
    nGA = ind_atual[35]
    nHF = ind_atual[36]
    nIG = ind_atual[37]
    nIH = ind_atual[38]
    nJI = ind_atual[39]

    solution = odeint(twoBody, Y0, dobra_pontos, args=(
    tauA, kAJ, int(nAJ), tauB, kBE, int(nBE), tauC, kCB, int(nCB), kCF, int(nCF), kCA, int(nCA), tauD, kDF, int(nDF),
    tauE, kEJ, int(nEJ), tauF, kFA, int(nFA), tauG, kGB, int(nGB), kGF, int(nGF), kGA, int(nGA), tauH, kHF, int(nHF),
    tauI, kIG, int(nIG), kIH, int(nIH), tauJ, kJI, int(nJI),))
    pA, pB, pC, pD, pE, pF, pG, pH, pI, pJ = organiza_pontos(solution)
    print(" ".join([str(x) for x in pA]))
    print(" ".join([str(x) for x in pB]))
    print(" ".join([str(x) for x in pC]))
    print(" ".join([str(x) for x in pD]))
    print(" ".join([str(x) for x in pE]))
    print(" ".join([str(x) for x in pF]))
    print(" ".join([str(x) for x in pG]))
    print(" ".join([str(x) for x in pH]))
    print(" ".join([str(x) for x in pI]))
    print(" ".join([str(x) for x in pJ]))
    return calcula_diferenca(pA, pB, pC, pD, pE, pF, pG, pH, pI, pJ)


# tauA = 1.73
# tauB = 2
# tauC = 0.81
# tauD = 0.11
# tauE = 1.23
# tauF = 1.78
# tauG = 1.14
# tauH = 1.04
# tauI = 3.47
# tauJ = 3.21
# kAJ = 0.45 *
# kBE = 0.56 *
# kCB = 0.99
# kCF = 0.77
# kCA = 0.71
# kDF = 0.66 *
# kEJ = 0.46 *
# kFA = 0.48 *
# kGB = 0.66
# kGF = 0.99
# kGA = 0.85
# kHF = 0.61 *
# kIG = 0.55
# kIH = 0.46
# kJI = 0.17 *
# nAJ = 20 *
# nBE = 9 *
# nCB = 24
# nCF = 12
# nCA = 2
# nDF = 2 *
# nEJ = 6 *
# nFA = 4 *
# nGB = 7
# nGF = 24
# nGA = 2
# nHF = 7 *
# nIG = 21
# nIH = 20
# nJI = 3 *