import matplotlib.pyplot as plt
import numpy as np


# ------------- CHANGER LES PARAMÈTRES AVANT DE LANCER LE CODE -------------------- # 
#initialiser la géométrie
a, b, c, d, e, f = 3, 2, 4, 2, 0.2, 6
N = 5
dx = 0.1 # step


largeur = (2*a + (N+1)*(c/2) + (N-1)*(d/2)) # trouvé en analysant l'image 
Lx,Ly = largeur, f # grandeur de la grille en x et en y sans espacement autour du pm
nx, ny = int(Lx/dx), int(Ly/dx) # nombre de points en x et en y

x = np.linspace(0, Lx, nx)
y = np.linspace(0, Ly, ny)
X, Y = np.meshgrid(x, y)

#initialiser le potentiel à 0 partout
V = np.zeros((ny,nx)) # potentiel initial
bloqué = np.zeros((ny,nx), dtype=bool)


# rentrer les dynodes pour laisser un potentiel qui ne varie pas
# Cette fonction est en [y,x]
#f, b et e sont en x
# a,c,d sont en y
def placer_dynodes_bas(V, bloqué):
    for i in range(N//2 + N%2):
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
        bloqué[iy_start:iy_end, ix_start:ix_end] = True
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

# Relaxation aka ce qui n'existe pas en gph sauf si t'es camille et tu bois en criss
# La variation minimale est établie à 10^-5 parce que 
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

# calculer le gradient en x et y grâce à numpy
Ey, Ex = np.gradient(-V, dx, dx)

E_norm = np.sqrt(Ex**2 + Ey**2) #norme


#implémenter une fonction de déterminer x(t) d'un électron
#Conditions initiales applicables: x(t=0) et v(t=0)

#constantes physiques:
q = -1.602*e-19 # charge de l'électron
m = 9.109*e-31 # masse
# F(x,y,z) = qE(x,y,z)

# Il faut approximer le champ entre les cases existantes avec euler pour des valeurs comme 0.25 ou 4.287
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

def position_el(x0, y0, vx0, vy0, Ex, Ey, dx, dt, steps):
    x = [x0]
    y = [y0]





#afficher le champ produit
plt.contourf(X, Y, E_norm, levels=100, cmap='plasma')
plt.colorbar(label="|E| (V/m)")

saut = 2
#vecteur
plt.quiver(X[::saut, ::saut], Y[::saut, ::saut],
           Ex[::saut, ::saut], Ey[::saut, ::saut],
           color='white', scale=60000)

plt.title("Champ électrique dans le tube photomultiplicateur")
plt.xlabel("x (mm)")
plt.ylabel("y (mm)")
plt.axis('equal')
plt.tight_layout()
plt.savefig("champ_PM.png", dpi=300)
plt.show()
