

import numpy as np
import matplotlib.pyplot as plt


# ------------- CHANGER LES PARAMÈTRES AVANT DE LANCER LE CODE -------------------- #

def initialiser_geo():

    #initialiser la géométrie
    a, b, c, d, e, f = 5, 3, 6, 3, 0.2, 10
    n = 4
    dx = 0.1 # step

    largeur = (2*a + (n+1)*(c/2) + (n-1)*(d/2)) # trouvé en analysant l'image
    lx,ly = largeur, f # grandeur de la grille en x et en y sans espacement autour du pm
    nx, ny = int(lx/dx), int(ly/dx) # nombre de points en x et en y

    x = np.linspace(0, lx, nx)
    y = np.linspace(-ly/2, ly/2, ny)
    x_grid, y_grid = np.meshgrid(x, y)

    #initialiser le potentiel à 0 partout
    v = np.zeros((ny,nx)) # potentiel initial #-----------------------------
    bloquer = np.zeros((ny,nx), dtype=bool)

    return [a,b,c,d,e,f,n,dx,lx,ly,nx,ny,x,y,x_grid,y_grid,v,bloquer]
    # 0=a, 1=b, 2=c, 3=d, 4=e, 5=f, 6=N, 7=dx, 8=Lx, 9=Ly, 10=nx,
    # 11=ny, 12=x, 13=y, 14=x_grid, 15=y_grid, 16=V, 17=bloqué

# rentrer les dynodes pour laisser un potentiel qui ne varie pas
# Cette fonction est en [y,x]
# f, b et e sont en x
# a,c,d sont en y

def placer_dynodes_bas(v, bloquer, parametres):

    a, b, c, d, e = parametres[0:5]
    n, dx = parametres[6:8]

    for i in range(n//2 + n%2): # 11//2 = 5 + 1
        pot_dyn = (2*i+1) * 100
        vert_start = b #coordonnée verticale du début de la dynode
        vert_end = b + e # la fin
        horiz_start = a + i*(c+d)   #coordonnée horizontale du début de la dynode
        horiz_end = horiz_start + c # fin

        ix_start = int(horiz_start/dx) #donne l'indice de colonne correspondant sur la grille
        ix_end = int(horiz_end/dx) # même chose mais pour la fin
        iy_start = int(vert_start/dx) # début y
        iy_end = int(vert_end/dx) # fin y

        v[iy_start:iy_end, ix_start:ix_end] = pot_dyn
        bloquer[iy_start:iy_end, ix_start:ix_end] = True #-------------
    return v, bloquer

def placer_dynodes_haut(v, bloquer, parametres):

    a, b, c, d, e, f, n, dx = parametres[0:8]

    for i in range(n//2):

        pot_dyn = (2*(i+1)) * 100
        vert_start = f-b-e #coordonnée verticale du début de la dynode
        vert_end = vert_start + e # la fin
        horiz_start = a + (i+1)*c +d/2 + i*d - c/2 #coordonnée horizontale du début de la dynode
        horiz_end = horiz_start + c # fin

        ix_start = int(horiz_start/dx) #donne l'indice de colonne correspondant sur la grille
        ix_end = int(horiz_end/dx) # même chose mais pour la fin
        iy_start = int(vert_start/dx) # début y
        iy_end = int(vert_end/dx) # fin y

        v[iy_start:iy_end, ix_start:ix_end] = pot_dyn
        bloquer[iy_start:iy_end, ix_start:ix_end] = True

    return v, bloquer

# Relaxation
# La variation minimale est établie à 10^-5 parce que c'est souvent ça dans d'autres problèmes

def relaxation(v, bloquer, parametres, variation=1e-5, max_iter=1000000):

    nx, ny = parametres[10:12]

    for iteration in range(max_iter):

        v_old = v.copy() # pour la comparaison pour la tolérance

        for i in range(1, ny-1): #tous les points sauf les bords

            for j in range(1, nx-1): # meme chose en x

                if not bloquer[i, j]:

                    v[i, j] = 0.25 * (v_old[i+1, j] + v_old[i-1, j] +
                                     v_old[i, j+1] + v_old[i, j-1])

        chang_pot_max = np.max(np.abs(v - v_old))

        if chang_pot_max <= variation:

            print(f"Convergence atteinte après {iteration} itérations "\
                  "avec une tolérance de variation de {variation}")
            print("Nous avons atteint un potentiel qui ne varie presque plus " \
            "et le problème est considéré résolu")

            break

        if iteration == max_iter - 1:

            print("Attention !!!!!!!!!! ")
            print("Le maximum d'itérations a été atteint sans stabilisation, " \
            "donc le programme a été arrêté")

    return v, bloquer

if __name__ == "__main__":
    # Appel de la fonction d'affichage avec les paramètres initiaux

    PARAMETRES = initialiser_geo()

    X_GRID, Y_GRID, V, BLOQUER = PARAMETRES[14:18]

    V, BLOQUER = placer_dynodes_haut(V, BLOQUER, PARAMETRES)
    V, BLOQUER = placer_dynodes_bas(V, BLOQUER, PARAMETRES)

    #Afficher le PM
    RES, BLOQUER = relaxation(V, BLOQUER, PARAMETRES, variation=1e-3, max_iter=4000)
    CP = plt.contourf(X_GRID, Y_GRID, RES, levels=100, cmap="plasma")
    # On ajoute le contour des dynodes pour être propre

    plt.contour(X_GRID, Y_GRID, BLOQUER, levels=[0.5], colors='black', linewidths=1)
    plt.colorbar(CP, label="Potentiel (V)")
    plt.title("Potentiel électrique dans le tube PM")
    plt.xlabel("x (mm)")
    plt.ylabel("y (mm)")
    plt.axis('equal')
    plt.tight_layout()
    plt.savefig("Question1_potentiel.png", dpi=300)
    plt.show()
