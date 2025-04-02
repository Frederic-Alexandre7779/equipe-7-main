import matplotlib.pyplot as plt
import numpy as np



# ------------- CHANGER LES PARAMÈTRES AVANT DE LANCER LE CODE -------------------- # 
#initialiser la géométrie
a, b, c, d, e, f = 0.2, 1, 4, 2, 0.2, 6
N = 6
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
def relaxation(V, bloqué, variation=1e-5, max_iter=1000000):
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
    if 0 <= i < ny and 0 <= j < nx: # vérifier que c'est dans la grille que j'ai créée sinon ça marche pas
        return Ex[i, j], Ey[i, j]
    else:
        return 0, 0 # si on n'est pas dans la grille, il n'y a pas de champ

# -------------question 3 c) (On modifie position_el avec une détection de contact avec contact_dynode)-------------#

# ------------POsitionner les dynodes du haut et celles du bas avec leur valeur de potentiel ------#

def position_dynodes_bas(dx):
    dynodes_bas = []
    for i in range(N//2 + N%2):
        vert_start = b #coordonnée verticale du début de la dynode
        vert_end = b + e # la fin
        horiz_start = a + i*(c+d)   #coordonnée horizontale du début de la dynode
        horiz_end = horiz_start + c # fin
        dynodes_bas.append((vert_start, vert_end, horiz_start, horiz_end, (2*i+1)*100))  # dynodes bas
    return dynodes_bas


def position_dynodes_haut(dx):
    dynodes_haut = []
    for i in range(N//2):
        vert_start = f-b #coordonnée verticale du début de la dynode
        vert_end = vert_start + e # la fin
        horiz_start = a + (i+1)*c +d/2 + i*d - c/2 #coordonnée horizontale du début de la dynode
        horiz_end = horiz_start + c # fin
        dynodes_haut.append((vert_start, vert_end, horiz_start, horiz_end, (2*(i+1)) * 100)) # dynodes haut
    return dynodes_haut

#-----------------après ces 2 fonctions, on a 4 coordonnées des côtés des rectangles créé par les dynodes soit du haut ou du bas -----#
#-----------------maintenant on défini s'il y a un contact avec les dynodes ------------------#

def contact_dyn_bas(x, y, dynodes_bas):
    for (x0, x1, y0, y1, _) in dynodes_bas:
        if x0 <= x <= x1 and y0 <= y <= y1:
            return True
    return False

def contact_dyn_haut(x, y, dynodes_haut):
    for (x0, x1, y0, y1, _) in dynodes_haut:
        if x0 <= x <= x1 and y0 <= y <= y1:
            return True
    return False

# -------- là j'ai fait 2 fonctions qui dise s'il y a un contact avec un dynodes --------# 


def position_el(x0, y0, vx0, vy0, Ex, Ey, dx, dt, it_max, dynodes_bas, dynodes_haut):
    q = -1.602*e-19 # charge de l'électron
    m = 9.109*e-31 # masse
    #placer l'électron à t=0
    x = [x0] # position à t=0
    y = [y0]
    vx = vx0 # vitesse à t=0
    vy = vy0

    for _ in range(it_max):
        # méthode d'Euler du champ au point courant
        Ex_val, Ey_val = Eulerer_lechamp(x[-1], y[-1], Ex, Ey, dx)
        #littéralement Euler:
        ax = (q*Ex_val / m) # accélération en x
        ay = (q*Ey_val / m)  # même chose en y

        vx += ax * dt # changement infinitésimal de la vitesse
        vy += ay * dt # même chose en y
        
        y_new = float(y[-1] + vy * dt)
        x_new = float(x[-1] + vx * dt) # changement infinitésimal de la position

        # on regarde s'il y a un contact avec celles du bas:
        if contact_dyn_bas(x_new, y_new, dynodes_bas):
            y_new += 2 # rebondi de 2mm vers le haut
            vy = -vy # la vitesse part dans l'autre sens pour les prochains calculs
            print(f"Rebond dynode BAS à [x,y] = [{x_new}{y_new}]")
            print("La vitesse a été inversé et la position a AUGMENTÉ de 2mm")

        # on regarde s'il y a un contact avec celles du haut :
        if contact_dyn_haut(x_new, y_new, dynodes_haut):
            y_new -= 2
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


    return np.array(x), np.array(y)


# ---------------------------------------Question 3b/c --------------------------------------------#

x0 = 0.5
y0 = 0
vx0 = 2
vy0 = 0
dt = 0.0001
it_max = 10000
dynodes_bas = position_dynodes_bas(dx)
dynodes_haut = position_dynodes_haut(dx)
traj_x, traj_y = position_el(x0, y0, vx0, vy0, Ex, Ey, dx, dt, it_max, dynodes_bas, dynodes_haut)

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
