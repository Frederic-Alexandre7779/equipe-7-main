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
y = np.linspace(-Ly/2, Ly/2, ny)
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

# Relaxation aka ce qui n'existe pas en gph sauf si t'es camille et tu bois en criss
# La variation minimale est établie à 10^-5 parce que 
def relaxation(V, bloqué, variation=1e-5, max_iter=10000):
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
print(V)
Ey, Ex = np.gradient(-V, dx, dx)

E_norm = np.sqrt(Ex**2 + Ey**2) #norme

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