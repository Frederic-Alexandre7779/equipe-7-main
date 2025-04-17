

import numpy as np
import matplotlib.pyplot as plt
import question1_clean


# ------------- CHANGER LES PARAMÈTRES AVANT DE LANCER LE CODE -------------------- #

def champ_elec(res, parametres):
    # Faire les calculs pour la question 2 seulement

    dx = parametres[7]

    ey, ex = np.gradient(-res, dx, dx)
    e_norm = np.sqrt(ex**2 + ey**2) #norme

    return [ex, ey, e_norm]


if __name__ == "__main__":
    # Appel de la fonction d'affichage pour le numéro 2 seulement

    PARAMETRES = question1_clean.initialiser_geo()

    X_GRID, Y_GRID, V, BLOQUER = PARAMETRES[14:18]

    V, BLOQUER = question1_clean.placer_dynodes_haut(V, BLOQUER, PARAMETRES)
    V, BLOQUER = question1_clean.placer_dynodes_bas(V, BLOQUER, PARAMETRES)

    RES, BLOQUER = question1_clean.relaxation(V, BLOQUER, PARAMETRES, variation=1e-3, max_iter=4000)

    EX, EY, E_NORM = champ_elec(RES, PARAMETRES)

    #afficher le champ produit
    plt.figure(figsize=(10, 5))
    plt.title("Champ électrique dans le tube photomultiplicateur")
    plt.contourf(X_GRID, Y_GRID, E_NORM, levels=100, cmap='plasma')
    plt.colorbar(label="|E| (V/m)")
    SAUT = 2
    # vecteur
    plt.quiver(X_GRID[::SAUT, ::SAUT], Y_GRID[::SAUT, ::SAUT],
               EX[::SAUT, ::SAUT], EY[::SAUT, ::SAUT],
               color='white', scale=10000)
    plt.xlabel("x (mm)")
    plt.ylabel("y (mm)")
    plt.axis('equal')
    plt.tight_layout()
    plt.savefig("Question2_champelec.png", dpi=300)
    plt.show()
