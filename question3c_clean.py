
import numpy as np
import matplotlib.pyplot as plt
import question1_clean
import question2clean
import question3a_clean

#----question 3 c) (On modifie position_el avec une détection de contact avec contact_dynode)-----#

# ------------POsitionner les dynodes du haut et celles du bas avec leur valeur de potentiel ------#

def position_dynodes_bas(i, parametres):

    a, b, c, d, e, f = parametres[0:6]

    vert_start = -f*0.5 + b #coordonnée verticale du début de la dynode
    vert_end = vert_start + e # la fin
    horiz_start = a + i*(c + d)   #coordonnée horizontale du début de la dynode
    horiz_end = horiz_start + c # fin
    pot = (2*i+1) * 100

    #dynodes_bas = [vert_start, vert_end, horiz_start, horiz_end, pot]
    return [vert_start, vert_end, horiz_start, horiz_end, pot]


def position_dynodes_haut(i, parametres):

    a, b, c, d, e, f = parametres[0:6]

    vert_start = f*0.5 - b #coordonnée verticale du début de la dynode
    vert_end = vert_start - e # la fin
    horiz_start = a + (i+1)*c +d/2 + i*d - c/2 #coordonnée horizontale du début de la dynode
    horiz_end = horiz_start + c # fin
    pot = (2*(i+1)) * 100

    #dynodes_haut = [vert_start, vert_end, horiz_start, horiz_end, pot]
    return [vert_start, vert_end, horiz_start, horiz_end, pot]


# Après ces 2 fonctions, on a 4 coordonnées des côtés des rectangles -------------#
# créé par les dynodes soit du haut ou du bas ------------------------------------#


#-----------------maintenant on défini s'il y a un contact avec les dynodes ------------------#


def contact_dyn_bas(x_new, y_new, x_old, y_old, parametres):
    # Juste pour être sûr d'avoir les bonnes dimensions

    f, n = parametres[5:7]
    lx = parametres[8]

    if -f/2 < y_new < 0 and 0 < x_new < lx:
        # Regarde si les données envoyées par la fonction de l'électron sont dans le PM

        pente = (y_new - y_old) / (x_new - x_old) # pente de la droite entre les 2 points
        y_initiale = -(pente*x_old - y_old) # trouve le b de y = mx + b

        for i in range (n//2 + n%2):
            # On initialise les dynodes avec leurs positions en mm
            dynodes_bas = position_dynodes_bas(i, parametres)

            if y_new <= dynodes_bas[1]: # On s'assure que le nouveau point est en bas de la dynode

                if x_new <= dynodes_bas[3] and x_old >= dynodes_bas[2]:
                    # Regarde si le point en x est au milieu d'une dynode

                    x_dyn = (dynodes_bas[1] - y_initiale)/pente # Trouve le x en haut de la dynode
                    y_dyn = dynodes_bas[1] # on remet le y à la hauteur centrale de la dynode

                    return [True, x_dyn, y_dyn] # Retourne la nouvelle position et True

                if x_new >= dynodes_bas[2] and x_old <= dynodes_bas[2]:
                    # Regarde si la pente passe à gauche de la dynode

                    if (dynodes_bas[0] - y_initiale)/pente >= dynodes_bas[2]:
                        # Regarde si la dynode est croisée tout court

                        if (dynodes_bas[1] - y_initiale)/pente <= dynodes_bas[2]:
                            # Si la dynode est croisée en haut

                            x_dyn = dynodes_bas[2] # Trouve le x en haut de la dynode côté gauche
                            y_dyn = dynodes_bas[1] # on remet le y à la hauteur de la dynode

                            return [True, x_dyn, y_dyn] # Retourne la nouvelle position et True

                        if (dynodes_bas[1] - y_initiale)/pente >= dynodes_bas[2]:
                            # Si la dynode est croisée en bas

                            x_dyn = (dynodes_bas[1] - y_initiale)/pente # Trouve le x de la dynode
                            y_dyn = dynodes_bas[1] # on remet le y à la hauteur de la dynode

                            return [True, x_dyn, y_dyn] # Retourne la nouvelle position et True

                if x_new >= dynodes_bas[3] and x_old <= dynodes_bas[3]:
                    # Regarde si la pente passe à droite de la dynode

                    if (dynodes_bas[1] - y_initiale)/pente <= dynodes_bas[3]:
                        # Si la dynode est croisée en haut

                        x_dyn = (dynodes_bas[1] - y_initiale)/pente # Trouve le x de la dynode
                        y_dyn = dynodes_bas[1] # on remet le y à la hauteur de la dynode

                        return  [True, x_dyn, y_dyn] # Retourne la nouvelle position et True

    return [False, x_new, y_new] # Si il n'y a pas de contact, on retourne la position d'origine


def contact_dyn_haut(x_new, y_new, x_old, y_old, parametres):
    # Juste pour être sûr d'avoir les bonnes dimensions

    f, n = parametres[5:7]
    lx = parametres[8]

    if f/2 > y_new > 0 and 0 < x_new < lx:
        # Regarde si les données envoyées par la fonction de l'électron sont dans le PM

        pente = (y_new - y_old) / (x_new - x_old) # pente de la droite entre les 2 points
        y_initiale = -(pente*x_old - y_old) # trouve le b de y = mx + b

        for i in range (n//2):
            # On initialise les dynodes avec leurs positions en mm
            dynodes_haut = position_dynodes_haut(i, parametres)

            if y_new >= dynodes_haut[1]: # On s'assure que le nouveau point est en haut de la dynode

                if x_new <= dynodes_haut[3] and x_old >= dynodes_haut[2]:
                    # Regarde si le point en x est au milieu d'une dynode

                    x_dyn = (dynodes_haut[1] - y_initiale)/pente # Trouve le x de la dynode
                    y_dyn = dynodes_haut[1] # on remet le y à la hauteur centrale de la dynode

                    return [True, x_dyn, y_dyn] # Retourne la nouvelle position et True

                if x_new >= dynodes_haut[2] and x_old <= dynodes_haut[2]:
                 # Regarde si la pente passe à gauche de la dynode

                    if (dynodes_haut[0] - y_initiale)/pente >= dynodes_haut[2]:
                        # Regarde si la dynode est croisée tout court

                        if (dynodes_haut[1] - y_initiale)/pente >= dynodes_haut[2]:
                            # Si la dynode est croisée en haut

                            x_dyn = (dynodes_haut[1] - y_initiale)/pente # Trouve le x de la dynode
                            y_dyn = dynodes_haut[1] # on remet le y à la hauteur de la dynode

                            return [True, x_dyn, y_dyn] # Retourne la nouvelle position et True

                        if (dynodes_haut[1] - y_initiale)/pente <= dynodes_haut[2]:
                            # Si la dynode est croisée en bas

                            x_dyn = dynodes_haut[2] # Trouve le x en bas de la dynode côté gauche
                            y_dyn = dynodes_haut[1] # on remet le y à la hauteur de la dynode

                            return [True, x_dyn, y_dyn] # Retourne la nouvelle position et True

                if x_new >= dynodes_haut[3] and x_old <= dynodes_haut[3]:
                    # Regarde si la pente passe à droite de la dynode

                    if (dynodes_haut[1] - y_initiale)/pente <= dynodes_haut[3]:
                        # Si la dynode est croisée en bas

                        x_dyn = (dynodes_haut[1] - y_initiale)/pente # Trouve le x de la dynode
                        y_dyn = dynodes_haut[1] # on remet le y à la hauteur centrale de la dynode

                        return [True, x_dyn, y_dyn] # Retourne la nouvelle position et True

    return [False, x_new, y_new] # Si il n'y a pas de contact, on retourne la position d'origine

# -------- là j'ai fait 2 fonctions qui dise s'il y a un contact avec un dynodes --------#

def position_el(x0, y0, vx0, vy0, ex, ey, dt, it_max, parametres):#

    dx, lx, ly = parametres[7:10]

    q = -1.602*10e-19 # charge de l'électron
    m = 9.109*10e-31 # masse de l'électron

    #placer l'électron à t=0
    x = [x0] # position à t=0
    y = [y0]
    vx = vx0 # vitesse à t=0
    vy = vy0

    # méthode d'Euler du champ au point courant
    ex_val, ey_val = question3a_clean.eulerer_lechamp(x[0], y[0], ex, ey, dx)

    #littéralement Euler:
    ax = (q*ex_val) / m # accélération en x
    ay = (q*ey_val) / m # même chose en y

    vx += ax * dt # changement infinitésimal de la vitesse
    vy += ay * dt # même chose en y

    y_new = float(y0 + (vy * dt))
    x_new = float(x0 + (vx * dt)) # changement infinitésimal de la position

    x.append(x_new) # ajoute le x et y après rebond
    y.append(y_new)

    for _ in range(it_max):

        contact_bas = contact_dyn_bas(x_new, y_new, x[_], y[_], parametres)
        contact_haut = contact_dyn_haut(x_new, y_new, x[_], y[_], parametres)

        # on regarde s'il y a un contact avec celles du bas:
        if contact_bas[0] is True:

            x_new = contact_bas[1]
            y_new = contact_bas[2]

            x.append(x_new) # on ajoute le x et y au contact
            y.append(y_new)

            print(f"Rebond dynode BAS à [x,y] = [{x_new}, {y_new}]")
            print("La vitesse a été inversé et la position a AUGMENTÉ de 2mm")

            y_new += 2 # on rebondi de 2mm vers le haut
            vy = -vy # la vitesse part dans l'autre sens pour les prochains calculs

            x.append(x_new) # on ajoute le x et y au contact
            y.append(y_new)

        # on regarde s'il y a un contact avec celles du haut :
        if contact_haut[0] is True:

            x_new = contact_haut[1]
            y_new = contact_haut[2]

            x.append(x_new) # on ajoute le x et y au contact
            y.append(y_new)

            print(f"Rebond dynode HAUT à [x,y] = [{x_new}, {y_new}]")
            print("La vitesse a été inversé et la position a DIMINUÉ de 2mm")

            y_new -= 2 # on rebondi de 2mm vers le bas
            vy = -vy

            x.append(x_new) # on ajoute le x et y au contact
            y.append(y_new)

        ex_val, ey_val = question3a_clean.eulerer_lechamp(x_new, y_new, ex, ey, dx)

        ax = (q*ex_val) / m # accélération en x
        ay = (q*ey_val) / m # même chose en y

        vx += ax * dt # changement infinitésimal de la vitesse
        vy += ay * dt # même chose en y

        y_new = float(y_new + (vy * dt))
        x_new = float(x_new + (vx * dt)) # changement infinitésimal de la position

        x.append(x_new)
        y.append(y_new)

        if x_new > lx and abs(y_new) < ly/2: #Regarde si l'électron à finit son trajet
            print("L'électron à passé le photomultiplicateur au complet :)")
            break

        #Là je suis tanné de faire des test avec l'électron qui sacre son camp en dehors du PM
        #faiq je rajoute une condition qui l'oblige à rester dans le PM
        if not (0 <= x_new < lx and abs(y_new) < ly/2):
            print("L'électron a crissé son camp")
            break #enfin ça va être moins chiant

    return np.array(x), np.array(y)

if __name__ == "__main__":

    PARAMETRES = question1_clean.initialiser_geo()

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

    traj_x, traj_y = position_el(X0, Y0, VX0, VY0, EX, EY, DT, IT_MAX, PARAMETRES)


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
    plt.quiver(X_GRID[::SAUT, ::SAUT], Y_GRID[::SAUT, ::SAUT],
            EX[::SAUT, ::SAUT], EY[::SAUT, ::SAUT],
            color='white', scale=6000)
    plt.plot(traj_x, traj_y, 'y-', label="Trajectoire")
    plt.plot(X0, Y0, 'go', label="Départ (x0,y0)")
    plt.xlabel("x (mm)")
    plt.ylabel("y (mm)")
    plt.title("Trajectoire")
    plt.legend()
    plt.axis("equal")
    plt.savefig("Q3c_trajectoire.png", dpi=300)
    plt.show()

    print("Position finale :", traj_x[-1], traj_y[-1])
    print("Déplacement total :", traj_x[-1] - traj_x[0], traj_y[-1] - traj_y[0])
    print("Champ au point de départ :", question3a_clean.eulerer_lechamp(X0, X0, EX, EY, DX))
