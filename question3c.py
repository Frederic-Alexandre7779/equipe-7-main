import matplotlib.pyplot as plt
import numpy as np



# ------------- CHANGER LES PARAMÈTRES AVANT DE LANCER LE CODE -------------------- # 
#initialiser la géométrie
a, b, c, d, e, f = 2, 2, 3, 2, 0.2, 10
N = 4
dx = 0.01 # step


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
        vert_start = f-b #coordonnée verticale du début de la dynode
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

# Relaxation
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
c = relaxation(V, bloqué, variation=1e-3, max_iter=3000)
print(c)

# --------------------------------------- Question 2 ---------------------------------#

# calculer le gradient en x et y grâce à numpy
Ey, Ex = np.gradient(-c, dx, dx)

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

def position_dynodes_bas(i, a=2, b=2, c=3, d=2, e=0.2, f=10):
    #dynodes_bas = []
    vert_start = -f*0.5 + b #coordonnée verticale du début de la dynode
    vert_end = vert_start + e # la fin
    horiz_start = a + i*(c + d)   #coordonnée horizontale du début de la dynode
    horiz_end = horiz_start + c # fin
    pot = (2*i+1) * 100
    #dynodes_bas = [vert_start, vert_end, horiz_start, horiz_end, pot]  # dynodes bas
    return [vert_start, vert_end, horiz_start, horiz_end, pot]


def position_dynodes_haut(i, a=2, b=2, c=3, d=2, e=0.2, f=10):
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

def contact_dyn_bas(x_new, y_new, x_old, y_old):#, dynodes_bas):
    a, b, c, d, e, f = 2, 2, 3 ,2 ,0.2, 10
    if y_new < 0 and 0 < x_new < (2*a + (N+1)*(c/2) + (N-1)*(d/2)):
        pente = (y_new - y_old) / (x_new - x_old)
        for i in range (N//2 + N%2):
            dynodes_bas = position_dynodes_bas(i, a, b, c, d, e, f)
            if x_new <= dynodes_bas[3] and x_old >= dynodes_bas[2]:
                if y_new <= dynodes_bas[1]:
                    x_dyn = dynodes_bas[1]/pente + x_new
                    y_dyn = dynodes_bas[1]
                    return [True, x_dyn, y_dyn]
            if x_new >= dynodes_bas[2] and x_old <= dynodes_bas[2]:
                if dynodes_bas[1]/pente >= dynodes_bas[2]:
                    x_dyn = dynodes_bas[1]/pente + x_new
                    y_dyn = dynodes_bas[1]
                    return [True, x_dyn, y_dyn]
            if x_new >= dynodes_bas[3] and x_old <= dynodes_bas[3]:
                if dynodes_bas[1]/pente <= dynodes_bas[3]:
                    x_dyn = dynodes_bas[1]/pente + x_new
                    y_dyn = dynodes_bas[1]
                    return  [True, x_dyn, y_dyn]
    return [False, x_new, y_new]

def contact_dyn_haut(x_new, y_new, x_old, y_old):#, dynodes_haut):
    a, b, c, d, e, f = 2, 2, 3 ,2 ,0.2, 10
    if y_new > 0 and 0 < x_new < (2*a + (N+1)*(c/2) + (N-1)*(d/2)):
        pente = y_new - y_old / x_new - x_old
        for i in range (N//2):
            dynodes_haut = position_dynodes_haut(i, a, b, c, d, e, f)
            if x_new <= dynodes_haut[3] and x_old >= dynodes_haut[2]:
                if y_new >= dynodes_haut[1]:
                    x_dyn = dynodes_haut[1]/pente + x_new
                    y_dyn = dynodes_haut[1]
                    return [True, x_dyn, y_dyn]
            if x_new >= dynodes_haut[2] and x_old <= dynodes_haut[2]:
                if dynodes_haut[1]/pente >= dynodes_haut[2]:
                    x_dyn = dynodes_haut[1]/pente + x_new
                    y_dyn = dynodes_haut[1]
                    return [True, x_dyn, y_dyn]
            if x_new >= dynodes_haut[3] and x_old <= dynodes_haut[3]:
                if dynodes_haut[1]/pente <= dynodes_haut[3]:
                    x_dyn = dynodes_haut[1]/pente + x_new
                    y_dyn = dynodes_haut[1]
                    return [True, x_dyn, y_dyn]
    return [False, x_new, y_new]

# -------- là j'ai fait 2 fonctions qui dise s'il y a un contact avec un dynodes --------# 


def position_el(x0, y0, vx0, vy0, Ex, Ey, dx, dt, it_max):#, dynodes_bas, dynodes_haut):
    q = -1.602*10e-19 # charge de l'électron
    m = 9.109*10e-31 # masse
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
        print(x_new, y_new)
        Ex_val, Ey_val = Eulerer_lechamp(x_new, y_new, Ex, Ey, dx) #change la valeur du champ mais ne marche pas en x

        ax = (q*Ex_val / m) # accélération en x
        ay = (q*Ey_val / m) # même chose en y

        vx += ax * dt # changement infinitésimal de la vitesse
        vy += ay * dt # même chose en y

        contact_bas = contact_dyn_bas(x_new, y_new, x[_], y[_])
        contact_haut = contact_dyn_haut(x_new, y_new, x[_], y[_])

        # on regarde s'il y a un contact avec celles du bas:
        if contact_bas[0] is True:#, dynodes_bas):
            x_new = contact_bas[1]
            y_new = contact_bas[2] + 2 # rebondi de 2mm vers le haut
            vy = -vy # la vitesse part dans l'autre sens pour les prochains calculs
            print(f"Rebond dynode BAS à [x,y] = [{x_new}{y_new}]")
            print("La vitesse a été inversé et la position a AUGMENTÉ de 2mm")

        # on regarde s'il y a un contact avec celles du haut :
        if contact_haut[0] is True:#, dynodes_haut):
            x_new = contact_haut[1]
            y_new = contact_haut[2] - 2
            vy = -vy
            print(f"Rebond dynode HAUT à [x,y] = [{x_new}{y_new}]")
            print("La vitesse a été inversé et la position a DIMINUÉ de 2mm")

        x.append(x_new)
        y.append(y_new)

        #Là je suis tanné de faire des test avec l'électron qui sacre son camp en dehors du PM
        #faiq je rajoute une condition qui l'oblige à rester dans le PM
        if not (0 <= x_new < Lx and -Ly/2 <= y_new <= Ly/2):
            print("L'électron a crissé son camp")
            break #enfin ça va être moins chiant

        y_new = float(y_new + (vy * dt))
        x_new = float(x_new + (vx * dt)) # changement infinitésimal de la position

    return np.array(x), np.array(y)


# ---------------------------------------Question 3b/c --------------------------------------------#

x0 = 2
y0 = 0
vx0 = 3
vy0 = 0
dt = 0.0001
it_max = 10000
#dynodes_bas = position_dynodes_bas()
#dynodes_haut = position_dynodes_haut()
traj_x, traj_y = position_el(x0, y0, vx0, vy0, Ex, Ey, dx, dt, it_max)#, dynodes_bas, dynodes_haut)

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
