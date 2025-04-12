import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Constantes
L = 1  # longueur du domaine
c = 1  # vitesse

# Fonction u(x,t)
def u(x, t, N_max):
    result = np.zeros_like(x)
    for n in range(1, N_max + 1):  # n entier seulement
        coeff = (8 * L) / (n**2 * np.pi**2)
        result += coeff * np.cos(n * np.pi * c * t / (2 * L)) * np.cos(n * np.pi * x / (2 * L))
    return result

# Grille de points
x = np.linspace(0, 10, 500)

# Initialisation du graphique
fig, ax = plt.subplots(figsize=(10, 6))
line, = ax.plot([], [], lw=2)
ax.set_xlim(0, 10)
ax.set_ylim(-1, 1)
ax.set_xlabel('x')
ax.set_ylabel('u(x,t)')
ax.set_title('Animation de u(x,t) avec N=2')
ax.grid(True)

# Fonction d'initialisation
def init():
    line.set_data([], [])
    return line,

plt.show()