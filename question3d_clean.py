
import numpy as np
import matplotlib.pyplot as plt
import question1_clean
import question2clean
import question3a_clean
import question3c_clean

#----question 3 d) (Fait le chemin de l'électron pour 12 dynodes)-----#

if __name__ == "__main__":

    A, B, C, D, E, F = 5, 3, 6, 3, 0.2, 10
    N = 12
    DX = 0.1 # step

    LARGEUR = (2*A + (N+1)*(C/2) + (N-1)*(D/2)) # trouvé en analysant l'image
    LX, LY = LARGEUR, F # grandeur de la grille en x et en y sans espacement autour du pm
    NX, NY = int(LX/DX), int(LY/DX) # nombre de points en x et en y

    X = np.linspace(0, LX, NX)
    Y = np.linspace(-LY/2, LY/2, NY)
    X_GRID, Y_GRID = np.meshgrid(X, Y)

    #initialiser le potentiel à 0 partout
    V = np.zeros((NY,NX)) # potentiel initial #-----------------------------
    BLOQUER = np.zeros((NY,NX), dtype=bool)

    PARAMETRES = A, B, C, D, E, F, N, DX, LX, LY, NX, NY, X, Y, X_GRID, Y_GRID, V, BLOQUER

    # On refait la géométrie parce qu'on veut 12 dynodes dans le problème.
    # Ça évite de changer le nombre de dynodes pour chaque question.

    DX, LX, LY = PARAMETRES[7:10]

    X_GRID, Y_GRID, V, BLOQUER = PARAMETRES[14:18]

    V, BLOQUER = question1_clean.placer_dynodes_haut(V, BLOQUER, PARAMETRES)
    V, BLOQUER = question1_clean.placer_dynodes_bas(V, BLOQUER, PARAMETRES)

    RES, BLOQUER = question1_clean.relaxation(V, BLOQUER, PARAMETRES, variation=1e-3, max_iter=4000)

    EX, EY, E_NORM = question2clean.champ_elec(RES, PARAMETRES)[0:5]

    X0 = 0
    Y0 = 0
    VX0 = 0
    VY0 = 0
    DT = 0.00000001
    IT_MAX = 50000

    traj_x, traj_y = question3c_clean.position_el(X0, Y0, VX0, VY0, EX, EY, DT, IT_MAX, PARAMETRES)


    fig, axe = plt.subplots(figsize=(10, 6))

    # Fond du champ électrique
    CP = axe.contourf(X_GRID, Y_GRID, E_NORM, levels=100, cmap='plasma')
    plt.colorbar(CP, label="|E| (V/m)")

    # Champ électrique (vecteurs blancs)
    SAUT = 2
    axe.quiver(X_GRID[::SAUT, ::SAUT], Y_GRID[::SAUT, ::SAUT], EX[::SAUT, ::SAUT],
               EY[::SAUT, ::SAUT], color='white', scale=6000)

    # Ligne de la trajectoire (vide au début)
    line, = axe.plot([], [], 'y-', label="Trajectoire de l'électron")
    point, = axe.plot([], [], 'go')  # point de l'électron

    # Mise en forme du graphique
    axe.set_xlim(0, LX)
    axe.set_ylim(-LY/2, LY/2)
    axe.set_xlabel("x (mm)")
    axe.set_ylabel("y (mm)")
    axe.set_title("Animation de la trajectoire de l'électron")
    axe.legend()
    axe.axis('equal')

    # Fonction pour mettre à jour l'animation à chaque frame
    def update(frame):
        line.set_data(traj_x[:frame], traj_y[:frame])  # Trace la trajectoire jusqu'au frame courant
        point.set_data(traj_x[frame-1], traj_y[frame-1])  # Place le point sur la position actuelle
        return line, point

    # Créer l'animation
    #ani = animation.FuncAnimation(fig, update, frames=len(traj_x), interval=10, blit=True)

    # Sauvegarder en mp4 (optionnel)
    # ani.save('animation_trajectoire.mp4', writer='ffmpeg', fps=60)

    #-------------------------------affichage de la trajectoire x(t)-------------------------------#

    #afficher le champ produit
    plt.contourf(X_GRID, Y_GRID, E_NORM, levels=100, cmap='plasma')
    plt.contour(X_GRID, Y_GRID, BLOQUER, levels=[0.5], colors='black', linewidths=1)

    SAUT = 2
    #vecteur
    #plt.quiver(X_GRID[::SAUT, ::SAUT], Y_GRID[::SAUT, ::SAUT],
            #EX[::SAUT, ::SAUT], EY[::SAUT, ::SAUT],
            #color='white', scale=6000)
    plt.plot(traj_x, traj_y, 'y-', label="Trajectoire")
    plt.plot(X0, Y0, 'go', label="Départ (x0,y0)")
    plt.xlabel("x (mm)")
    plt.ylabel("y (mm)")
    plt.title("Trajectoire")
    plt.legend()
    plt.axis("equal")
    plt.savefig("Q3d_trajectoire.png", dpi=300)
    plt.show()

    print("Position finale :", traj_x[-1], traj_y[-1])
    print("Déplacement total :", traj_x[-1] - traj_x[0], traj_y[-1] - traj_y[0])
    print("Champ au point de départ :", question3a_clean.eulerer_lechamp(X0, X0, EX, EY, DX))
