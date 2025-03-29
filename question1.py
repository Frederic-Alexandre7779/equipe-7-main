#question 1:
#méthode de la relaxation pour calculer le potentiel en tout point
#1ere à 0V, 2e à 100V etc
#Trouver un critère d'arrêt dans le calcul
import matplotlib as mpl
import math
import numpy as np

def construire_pm(N, a, b, c, d, e, f):
    dimensions_pm = (2*a + (N+1)c/2 + (N-1)d/2 ,f)
    pot_dyn = []
    for i in range(1, N+1):
        pot_dyn.append(i*100)

    res = []
    print(pot_dyn)
    return res
    



print(construire_pm(N=4, a=3, b=2, c=4, d=2, e=0.2, f=6))