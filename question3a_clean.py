

import numpy as np
import question1_clean
import question2clean

# ------------- CHANGER LES PARAMÈTRES AVANT DE LANCER LE CODE -------------------- #

# ---------------------------------------- Question 3a--------------------------------#

#implémenter une fonction de déterminer x(t) d'un électron
#Conditions initiales applicables: x(t=0) et v(t=0)

#fonctions importantes:
# F(x,y,z) = qE(x,y,z)
# a = qE/m

# Il faut approximer le champ entre les cases existantes
# avec euler pour des valeurs comme 0.25 ou 4.287

def eulerer_lechamp(x_p, y_p, ex, ey, dx):

    # transformer la position en indice de row/colonne
    i = int(y_p / dx)
    j = int(x_p / dx)

    # il faut créer une référence de grandeur pour savoir si on est dans la grille
    ny, nx = ex.shape

    ix = 0 #Détermine la position en x par rapport au subdivisions total ex: 35/60

    if i <= 0: #Pour isoler la hauteur de 0 à 60 au lieu de -30 à 30

        ix = int(i + ny*0.5)

    if i > 0:

        ix = int(i + ny*0.5)

    if 0 <= ix <= ny and 0 <= j <= nx:
        # vérifier que c'est dans la grille que j'ai créée sinon ça marche pas

        return ex[ix, j], ey[ix, j]

    else:

        return 0, 0 # si on n'est pas dans la grille, il n'y a pas de champ

def position_el(x_ini, y_ini, vx_ini, vy_ini, ex, ey, step, max_it, parametres):

    q = -1.602e-19 # charge de l'électron
    m = 9.109e-31 # masse

    #placer l'électron à t=0
    x = [x_ini] # position à t=0
    y = [y_ini]
    vx = vx_ini # vitesse à t=0
    vy = vy_ini

    dx = parametres[7]
    lx, ly = parametres[8:10]

    ex_val, ey_val = eulerer_lechamp(x_ini, y_ini, ex, ey, dx) # champ au point de départ

    ax = (q*ex_val) / m # accélération en x
    ay = (q*ey_val) / m  # même chose en y

    vx += ax * step # changement infinitésimal de la vitesse
    vy += ay * step # même chose en y

    x_new = float(x_ini + vx * step) # changement infinitésimal de la position
    y_new = float(y_ini + vy * step)

    x.append(x_new)
    y.append(y_new)

    for _ in range(max_it):

        # méthode d'Euler du champ au point courant
        ex_val, ey_val = eulerer_lechamp(x_new, y_new, ex, ey, dx)

        #littéralement Euler:
        ax = (q*ex_val) / m # accélération en x
        ay = (q*ey_val) / m  # même chose en y

        vx += ax * step # changement infinitésimal de la vitesse
        vy += ay * step # même chose en y

        x_new = float(x_new + vx * step) # changement infinitésimal de la position
        y_new = float(y_new + vy * step)

        x.append(x_new)
        y.append(y_new)

        if not (0 <= x_new < lx and abs(y_new) <= ly/2):

            print("L'électron a crissé son camp")
            break #enfin ça va être moins chiant

    print("Champ au point de départ :", eulerer_lechamp(x_ini, y_ini, ex, ey, dx))

    return np.array(x), np.array(y)

if __name__ == "__main__":
    # Appel de la fonction d'affichage pour le numéro 3a seulement

    PARAMETRES = question1_clean.initialiser_geo()

    X0 = 0 # position initiale en mm
    Y0 = 0 # position initiale en mm
    VX0 = 0 # vitesse initiale en mm/s
    VY0 = 0 # vitesse initiale en mm/s
    DT = 1e-9 # pas de temps en s
    IT_MAX = 10000 # nombre d'itérations maximum

    X_GRID, Y_GRID, V, BLOQUER = PARAMETRES[14:18]

    V, BLOQUER = question1_clean.placer_dynodes_haut(V, BLOQUER, PARAMETRES)
    V, BLOQUER = question1_clean.placer_dynodes_bas(V, BLOQUER, PARAMETRES)

    RES, BLOQUER = question1_clean.relaxation(V, BLOQUER, PARAMETRES, variation=1e-3, max_iter=3000)

    EX, EY = question2clean.champ_elec(RES, PARAMETRES)[0:2]

    traj_x, traj_y = position_el(X0, Y0, VY0, VY0, EX, EY, DT, IT_MAX, PARAMETRES)
