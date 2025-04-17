import numpy as np
import matplotlib.pyplot as plt
#from matplotlib import animation



# ------------- CHANGER LES PARAMÈTRES AVANT DE LANCER LE CODE -------------------- #
#initialiser la géométrie
a, b, c, d, e, f = 4, 3, 6, 4, 0.2, 10
N = 4
dx = 0.1 # step


largeur = (2*a + (N+1)*(c/2) + (N-1)*(d/2)) # trouvé en analysant l'image
Lx,Ly = largeur, f # grandeur de la grille en x et en y sans espacement autour du pm
nx, ny = int(Lx/dx), int(Ly/dx) # nombre de points en x et en y

x = np.linspace(0, Lx, nx)
y = np.linspace(-Ly/2, Ly/2, ny)
X, Y = np.meshgrid(x, y)

#initialiser le potentiel à 0 partout
V = np.zeros((ny,nx)) # potentiel initial, car c'est Laplace
bloqué = np.zeros((ny,nx), dtype=bool)


# rentrer les dynodes 
# Cette fonction est en [y,x]
#f, b et e sont en x
# a,c,d sont en y
def placer_dynodes_bas(V, bloqué):
    for i in range(N//2 + N%2): # 11//2 = 5 + 1
        pot_dyn = (2*i+1) * 100
        vert_start = b #coordonnée verticale du début de la dynode
        vert_end = b + e # la fin
        horiz_start = a + i*(c+d)   #coordonnée horizontale du début de la dynode
        horiz_end = horiz_start + c # fin

        ix_start = int(horiz_start/dx) #donne l'indice de colonne correspondant sur la grille
        ix_end = int(horiz_end/dx) # même chose mais pour la fin
        iy_start = int(vert_start/dx) # début y
        iy_end = int(vert_end/dx) # fin y

        V[iy_start:iy_end, ix_start:ix_end] = pot_dyn
        bloqué[iy_start:iy_end, ix_start:ix_end] = True #-------------
    return V, bloqué

def placer_dynodes_haut(V, bloqué):
    for i in range(N//2):
        pot_dyn = (2*(i+1)) * 100
        vert_start = f-b-e #coordonnée verticale du début de la dynode
        vert_end = vert_start + e # la fin
        horiz_start = a + (i+1)*c +d/2 + i*d - c/2 #coordonnée horizontale du début de la dynode
        horiz_end = horiz_start + c # fin

        ix_start = int(horiz_start/dx) #donne l'indice de colonne correspondant sur la grille
        ix_end = int(horiz_end/dx) # même chose mais pour la fin
        iy_start = int(vert_start/dx) # début y
        iy_end = int(vert_end/dx) # fin y

        V[iy_start:iy_end, ix_start:ix_end] = pot_dyn
        bloqué[iy_start:iy_end, ix_start:ix_end] = True
    return V, bloqué

V, bloqué = placer_dynodes_haut(V, bloqué)
V, bloqué = placer_dynodes_bas(V, bloqué)

# La variation minimale est établie à 10^-5 parce que c'est souvent ça dans d'autres problèmes
def relaxation(V, bloqué, variation= 1e-5, max_iter=1000000):
    for iteration in range(max_iter):
        V_old = V.copy() # pour la comparaison pour la tolérance
        for i in range(1, ny-1): #tous les points sauf les bords
            for j in range(1, nx-1): # meme chose en x
                if not bloqué[i, j]:
                    V[i, j] = 0.25 * (V_old[i+1, j] + V_old[i-1, j] +
                                     V_old[i, j+1] + V_old[i, j-1])
        chang_pot_max = np.max(np.abs(V - V_old))
        if chang_pot_max < variation:
            print(f"Convergence atteinte après {iteration} itérations avec une tolérance de variation de {variation})")
            print("Nous avons atteint un potentiel qui ne varie presque plus et le problème est considéré résolu")
            break
    else:
        print("Attention !!!!!!!!!! ")
        print("Le maximum d'itérations a été atteint sans stabilisation, donc le programme a été arrêté")
    return V

W = relaxation(V, bloqué, variation=1e-3, max_iter=3000)
print(W)

# --------------------------------------- Question 2 ---------------------------------#

# calculer le gradient en x et y grâce à numpy
Ey, Ex = np.gradient(-W, dx, dx)

E_norm = np.sqrt(Ex**2 + Ey**2) #norme

# ---------------------------------------- Question 3a--------------------------------#

# Il faut approximer le champ entre les cases existantes avec euler pour des valeurs 
def Eulerer_lechamp(x, y, Ex, Ey, dx):
    # transformer la position en indice de row/colonne
    i = int(y / dx)
    j = int(x / dx)
    # il faut créer une référence de grandeur pour savoir si on est dans la grille
    ny, nx = Ex.shape

    ix = 0 #Détermine la position en x par rapport au subdivisions total ex: 35/60
    if i <= 0: #Pour isoler la hauteur de 0 à 60 au lieu de -30 à 30
        ix = int(i + ny*0.5)
    if i > 0:
        ix = int(i + ny*0.5)
    
    if 0 <= ix <= ny and 0 <= j <= nx: # vérifier que c'est dans la grille que j'ai créée sinon ça marche pas
        return Ex[ix, j], Ey[ix, j]
    else:
        return 0, 0 # si on n'est pas dans la grille, il n'y a pas de champ

# -------------question 3 c) (On modifie position_el avec une détection de contact avec contact_dynode)-------------#

# ------------POsitionner les dynodes du haut et celles du bas avec leur valeur de potentiel ------#

def position_dynodes_bas(i):
    #dynodes_bas = []
    vert_start = -f*0.5 + b #coordonnée verticale du début de la dynode
    vert_end = vert_start + e # la fin
    horiz_start = a + i*(c + d)   #coordonnée horizontale du début de la dynode
    horiz_end = horiz_start + c # fin
    pot = (2*i+1) * 100
    #dynodes_bas = [vert_start, vert_end, horiz_start, horiz_end, pot]  # dynodes bas
    return [vert_start, vert_end, horiz_start, horiz_end, pot]


def position_dynodes_haut(i):
    #dynodes_haut = []
    vert_start = f*0.5 - b #coordonnée verticale du début de la dynode
    vert_end = vert_start - e # la fin
    horiz_start = a + (i+1)*c +d/2 + i*d - c/2 #coordonnée horizontale du début de la dynode
    horiz_end = horiz_start + c # fin
    pot = (2*(i+1)) * 100
    #dynodes_haut = [vert_start, vert_end, horiz_start, horiz_end, pot] # dynodes haut
    return [vert_start, vert_end, horiz_start, horiz_end, pot]

#-----------------après ces 2 fonctions, on a 4 coordonnées des côtés des rectangles créé par les dynodes soit du haut ou du bas -----#


#-----------------maintenant on défini s'il y a un contact avec les dynodes ------------------#


def contact_dyn_bas(x_new, y_new, x_old, y_old):
    # Juste pour être sûr d,avoir les bonnes dimensions

    if -f/2 < y_new < 0 and 0 < x_new < (2*a + (N+1)*(c/2) + (N-1)*(d/2)):
        # Regarde si les données envoyées par la fonction de l'électron sont dans le PM

        pente = (y_new - y_old) / (x_new - x_old) # pente de la droite entre les 2 points
        y_initiale = -(pente*x_old - y_old) # trouve le b de y = mx + b

        for i in range (N//2 + N%2):
            # On initialise les dynodes avec leurs positions en mm
            dynodes_bas = position_dynodes_bas(i)

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

                            x_dyn = dynodes_bas[2] # Trouve le x en haut de la dynode (Frappe à gauche)
                            y_dyn = dynodes_bas[1] # on remet le y à la hauteur centrale de la dynode

                            return [True, x_dyn, y_dyn] # Retourne la nouvelle position et True

                        if (dynodes_bas[1] - y_initiale)/pente >= dynodes_bas[2]:
                            # Si la dynode est croisée en bas

                            x_dyn = (dynodes_bas[1] - y_initiale)/pente # Trouve le x en haut de la dynode
                            y_dyn = dynodes_bas[1] # on remet le y à la hauteur centrale de la dynode

                            return [True, x_dyn, y_dyn] # Retourne la nouvelle position et True
                    
                if x_new >= dynodes_bas[3] and x_old <= dynodes_bas[3]:
                    # Regarde si la pente passe à droite de la dynode

                    if (dynodes_bas[1] - y_initiale)/pente <= dynodes_bas[3]:
                        # Si la dynode est croisée en haut

                        x_dyn = (dynodes_bas[1] - y_initiale)/pente # Trouve le x en haut de la dynode
                        y_dyn = dynodes_bas[1] # on remet le y à la hauteur centrale de la dynode

                        return  [True, x_dyn, y_dyn] # Retourne la nouvelle position et True
                
    return [False, x_new, y_new] # Si il n'y a pas de contact, on retourne la position d'origine et False


def contact_dyn_haut(x_new, y_new, x_old, y_old):
    # Juste pour être sûr d'avoir les bonnes dimensions

    if f/2 > y_new > 0 and 0 < x_new < (2*a + (N+1)*(c/2) + (N-1)*(d/2)):
        # Regarde si les données envoyées par la fonction de l'électron sont dans le PM

        pente = (y_new - y_old) / (x_new - x_old) # pente de la droite entre les 2 points
        y_initiale = -(pente*x_old - y_old) # trouve le b de y = mx + b

        for i in range (N//2):
            # On initialise les dynodes avec leurs positions en mm
            dynodes_haut = position_dynodes_haut(i)

            if y_new >= dynodes_haut[1]: # On s'assure que le nouveau point est en haut de la dynode

                if x_new <= dynodes_haut[3] and x_old >= dynodes_haut[2]:
                    # Regarde si le point en x est au milieu d'une dynode

                    x_dyn = (dynodes_haut[1] - y_initiale)/pente # Trouve le x de croisement en haut de la dynode
                    y_dyn = dynodes_haut[1] # on remet le y à la hauteur centrale de la dynode

                    return [True, x_dyn, y_dyn] # Retourne la nouvelle position et True
                
                if x_new >= dynodes_haut[2] and x_old <= dynodes_haut[2]:
                 # Regarde si la pente passe à gauche de la dynode

                    if (dynodes_haut[0] - y_initiale)/pente >= dynodes_haut[2]:
                        # Regarde si la dynode est croisée tout court

                        if (dynodes_haut[1] - y_initiale)/pente >= dynodes_haut[2]:
                            # Si la dynode est croisée en haut

                            x_dyn = (dynodes_haut[1] - y_initiale)/pente # Trouve le x en haut de la dynode
                            y_dyn = dynodes_haut[1] # on remet le y à la hauteur centrale de la dynode

                            return [True, x_dyn, y_dyn] # Retourne la nouvelle position et True
                    
                        if (dynodes_haut[1] - y_initiale)/pente <= dynodes_haut[2]:
                            # Si la dynode est croisée en bas

                            x_dyn = dynodes_haut[2] # Trouve le x en bas de la dynode (Frappe à gauche)
                            y_dyn = dynodes_haut[1] # on remet le y à la hauteur centrale de la dynode

                            return [True, x_dyn, y_dyn] # Retourne la nouvelle position et True
                    
                if x_new >= dynodes_haut[3] and x_old <= dynodes_haut[3]:
                    # Regarde si la pente passe à droite de la dynode

                    if (dynodes_haut[1] - y_initiale)/pente <= dynodes_haut[3]:
                        # Si la dynode est croisée en bas

                        x_dyn = (dynodes_haut[1] - y_initiale)/pente # Trouve le x en bas de la dynode
                        y_dyn = dynodes_haut[1] # on remet le y à la hauteur centrale de la dynode

                        return [True, x_dyn, y_dyn] # Retourne la nouvelle position et True
                    
    return [False, x_new, y_new] # Si il n'y a pas de contact, on retourne la position d'origine et False

# -------- là j'ai fait 2 fonctions qui dise s'il y a un contact avec un dynodes --------# 


def position_el(x0, y0, vx0, vy0, Ex, Ey, dx, dt, it_max):#, dynodes_bas, dynodes_haut):
    q = -1.602*10e-19 # charge de l'électron
    m = 9.109*10e-31 # masse de l'électron
    #placer l'électron à t=0
    x = [x0] # position à t=0
    y = [y0]
    vx = vx0 # vitesse à t=0
    vy = vy0

    # méthode d'Euler du champ au point courant
    Ex_val, Ey_val = Eulerer_lechamp(x[0], y[0], Ex, Ey, dx)
    #littéralement Euler:
    ax = (q*Ex_val / m) # accélération en x
    ay = (q*Ey_val / m) # même chose en y

    vx += ax * dt # changement infinitésimal de la vitesse
    vy += ay * dt # même chose en y

    y_new = float(y0 + (vy * dt))
    x_new = float(x0 + (vx * dt)) # changement infinitésimal de la position

    for _ in range(it_max):

        print(x_new, y_new) # pour voir si ça bouge comme il faut

        Ex_val, Ey_val = Eulerer_lechamp(x_new, y_new, Ex, Ey, dx) #change la valeur du champ mais ne marche pas en x

        ax = (q*Ex_val / m) # accélération en x
        ay = (q*Ey_val / m) # même chose en y

        vx += ax * dt # changement infinitésimal de la vitesse
        vy += ay * dt # même chose en y

        contact_bas = contact_dyn_bas(x_new, y_new, x[_], y[_])
        contact_haut = contact_dyn_haut(x_new, y_new, x[_], y[_])

        x.append(x_new)
        y.append(y_new)

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

            x.append(x_new) # ajoute le x et y après rebond
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

            x.append(x_new) # on ajoute le x et y au rebond
            y.append(y_new)

            #Là je suis tanné de faire des test avec l'électron qui sacre son camp en dehors du PM
            #faiq je rajoute une condition qui l'oblige à rester dans le PM
            if not (0 <= x_new < largeur and abs(y_new) < Ly/2):
                print("L'électron a crissé son camp")
                break #enfin ça va être moins chiant

        y_new = float(y_new + (vy * dt))
        x_new = float(x_new + (vx * dt)) # changement infinitésimal de la position

        if x_new > largeur and abs(y_new) < Ly/2: #Regarde si l'électron à finit son trajet
            print("L'électron à passé le photomultiplicateur au complet :)")
            break 

        #Là je suis tanné de faire des test avec l'électron qui sacre son camp en dehors du PM
        #faiq je rajoute une condition qui l'oblige à rester dans le PM
        if not (0 <= x_new < largeur and abs(y_new) < Ly/2):
            print("L'électron a crissé son camp")
            break #enfin ça va être moins chiant
        
    return np.array(x), np.array(y)


# ---------------------------------------Question 3b/c --------------------------------------------#

x0 = 0
y0 = 0
vx0 = 0
vy0 = 0
dt = 0.00000001
it_max = 50000
#dynodes_bas = position_dynodes_bas()
#dynodes_haut = position_dynodes_haut()
traj_x, traj_y = position_el(x0, y0, vx0, vy0, Ex, Ey, dx, dt, it_max)#, dynodes_bas, dynodes_haut)


fig, axe = plt.subplots(figsize=(10, 6))

# Fond du champ électrique
cp = axe.contourf(X, Y, E_norm, levels=100, cmap='plasma')
plt.colorbar(cp, label="|E| (V/m)")

# Champ électrique (vecteurs blancs)
saut = 2
axe.quiver(X[::saut, ::saut], Y[::saut, ::saut], Ex[::saut, ::saut], Ey[::saut, ::saut], color='white', scale=6000)

# Ligne de la trajectoire (vide au début)
line, = axe.plot([], [], 'y-', label="Trajectoire de l'électron")
point, = axe.plot([], [], 'go')  # point de l'électron

# Mise en forme du graphique
axe.set_xlim(0, Lx)
axe.set_ylim(-Ly/2, Ly/2)
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

plt.show()


#----------------------------------------affichage de la trajectoire x(t) premier attempt------------------------#

#afficher le champ produit
plt.contourf(X, Y, E_norm, levels=100, cmap='plasma')
plt.colorbar(label="|E| (V/m)")

saut = 2
#vecteur
plt.quiver(X[::saut, ::saut], Y[::saut, ::saut],
            Ex[::saut, ::saut], Ey[::saut, ::saut],
            color='white', scale=6000)
plt.plot(traj_x, traj_y, 'y-', label="Trajectoire")
plt.plot(x0, y0, 'go', label="Départ (x0,y0)")
plt.xlabel("x (mm)")
plt.ylabel("y (mm)")
plt.title("Trajectoire")
plt.legend()
plt.axis("equal")
plt.savefig("Q3c_trajectoire.png", dpi=300)
plt.show()

print("Position finale :", traj_x[-1], traj_y[-1])
print("Déplacement total :", traj_x[-1] - traj_x[0], traj_y[-1] - traj_y[0])
print("Champ au point de départ :", Eulerer_lechamp(x0, y0, Ex, Ey, dx))
