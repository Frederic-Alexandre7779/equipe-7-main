import numpy as np
import matplotlib.pyplot as plt

# Paramètres géométriques du tube PM
a, b, c, d, e, f = 3, 2, 4, 2, 0.2, 6  # en mm
N_dynodes = 4
V0 = 100  # Différence de potentiel entre dynodes
Lx, Ly = 12, 12  # Taille du domaine (mm)
dx = 0.1  # Pas spatial (mm)
nx, ny = int(Lx/dx), int(Ly/dx)

# Grille de potentiel
V = np.zeros((ny, nx))

# Coordonnées de la grille
x = np.linspace(0, Lx, nx)
y = np.linspace(0, Ly, ny)
X, Y = np.meshgrid(x, y)

# Fonction utilitaire pour placer les dynodes
def placer_dynodes(V):
    for i in range(N_dynodes):
        V_dyn = (i + 1) * V0
        y_pos = a + i * d
        x_start = (f - c) / 2
        x_end = x_start + c
        ix_start = int(x_start / dx)
        ix_end = int(x_end / dx)
        iy = int(y_pos / dx)
        V[iy, ix_start:ix_end] = V_dyn
    return V

# Appliquer les conditions aux limites
V = placer_dynodes(V)

# Masque pour les points fixes (dynodes et parois)
masque_fixe = (V != 0)

# Méthode de relaxation
def relaxation(V, tol=1e-5, max_iter=10000):
    for n in range(max_iter):
        V_old = V.copy()
        for i in range(1, ny-1):
            for j in range(1, nx-1):
                if not masque_fixe[i, j]:
                    V[i, j] = 0.25 * (V_old[i+1, j] + V_old[i-1, j] +
                                     V_old[i, j+1] + V_old[i, j-1])
        delta = np.max(np.abs(V - V_old))
        if delta < tol:
            print(f"Convergence atteinte après {n} itérations (Δ = {delta:.2e})")
            break
    else:
        print("Avertissement : nombre maximal d'itérations atteint.")
    return V

# Calcul du potentiel
V = relaxation(V)

# Affichage
plt.figure(figsize=(8, 6))
cp = plt.contourf(X, Y, V, levels=50, cmap='inferno')
plt.colorbar(cp, label='Potentiel (V)')
plt.title('Potentiel électrique dans le tube PM')
plt.xlabel('x (mm)')
plt.ylabel('y (mm)')
plt.axis('equal')
plt.tight_layout()
plt.savefig("potentiel_PM.png")
plt.show()
