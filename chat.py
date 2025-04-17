import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation

# --------- PARAMÈTRES DU PROBLÈME ---------
a, b, c, d, e, f = 0.2, 2, 4, 2, 0.2, 6
N = 4
dx = 0.1

largeur = (2*a + (N+1)*(c/2) + (N-1)*(d/2))
Lx = largeur
Ly = f
nx = int(Lx/dx)
ny = int(Ly/dx)

x = np.linspace(0, Lx, nx)
y = np.linspace(-Ly/2, Ly/2, ny)
X, Y = np.meshgrid(x, y)

V = np.zeros((ny, nx))
bloqué = np.zeros((ny, nx), dtype=bool)

# --------- DÉFINITION DES DYNODES ---------
def placer_dynodes_bas(V, bloqué):
    for i in range(N//2 + N%2):
        pot_dyn = (2*i+1) * 100
        vert_start = b
        vert_end = b + e
        horiz_start = a + i*(c+d)
        horiz_end = horiz_start + c

        ix_start = int(horiz_start/dx)
        ix_end = int(horiz_end/dx)
        iy_start = int(vert_start/dx)
        iy_end = int(vert_end/dx)

        V[iy_start:iy_end, ix_start:ix_end] = pot_dyn
        bloqué[iy_start:iy_end, ix_start:ix_end] = True
    return V, bloqué

def placer_dynodes_haut(V, bloqué):
    for i in range(N//2):
        pot_dyn = (2*(i+1)) * 100
        vert_start = f-b
        vert_end = vert_start + e
        horiz_start = a + (i+1)*c + d/2 + i*d - c/2
        horiz_end = horiz_start + c

        ix_start = int(horiz_start/dx)
        ix_end = int(horiz_end/dx)
        iy_start = int(vert_start/dx)
        iy_end = int(vert_end/dx)

        V[iy_start:iy_end, ix_start:ix_end] = pot_dyn
        bloqué[iy_start:iy_end, ix_start:ix_end] = True
    return V, bloqué

V, bloqué = placer_dynodes_haut(V, bloqué)
V, bloqué = placer_dynodes_bas(V, bloqué)

# --------- RELAXATION ---------
def relaxation(V, bloqué, variation=1e-5, max_iter=1000000):
    for iteration in range(max_iter):
        V_old = V.copy()
        for i in range(1, ny-1):
            for j in range(1, nx-1):
                if not bloqué[i, j]:
                    V[i, j] = 0.25 * (V_old[i+1, j] + V_old[i-1, j] + V_old[i, j+1] + V_old[i, j-1])
        chang_pot_max = np.max(np.abs(V - V_old))
        if chang_pot_max < variation:
            print(f"Convergence atteinte après {iteration} itérations")
            break
    else:
        print("Maximum d'itérations atteint sans convergence")
    return V

res = relaxation(V, bloqué)

Ey, Ex = np.gradient(-res, dx, dx)
E_norm = np.sqrt(Ex**2 + Ey**2)

# --------- OUTILS SEGMENT ---------
def ccw(A, B, C):
    return (C[1]-A[1]) * (B[0]-A[0]) > (B[1]-A[1]) * (C[0]-A[0])

def segments_intersect(A, B, C, D):
    return (ccw(A, C, D) != ccw(B, C, D)) and (ccw(A, B, C) != ccw(A, B, D))

def segment_intersect_rectangle(x_old, y_old, x_new, y_new, x_start, x_end, y_start, y_end):
    A = (x_old, y_old)
    B = (x_new, y_new)
    rect_sides = [
        ((x_start, y_start), (x_end, y_start)),
        ((x_end, y_start), (x_end, y_end)),
        ((x_end, y_end), (x_start, y_end)),
        ((x_start, y_end), (x_start, y_start))
    ]
    for C, D in rect_sides:
        if segments_intersect(A, B, C, D):
            return True
    return False

# --------- MOUVEMENT ÉLECTRON ---------
q = -1.602e-19
m = 9.109e-31

def Eulerer_lechamp(x, y, Ex, Ey, dx):
    i = int(y / dx)
    j = int(x / dx)
    ny, nx = Ex.shape
    if 0 <= i < ny and 0 <= j < nx:
        return Ex[i, j], Ey[i, j]
    else:
        return 0, 0

def position_el(x0, y0, vx0, vy0, Ex, Ey, dx, dt, it_max, a, b, c, d, e, f, N):
    x = [x0]
    y = [y0]
    vx = vx0
    vy = vy0

    positions_dynodes_bas = []
    positions_dynodes_haut = []

    for i in range(N//2 + N%2):
        x_start = a + i*(c+d)
        x_end = x_start + c
        y_start = b
        y_end = b + e
        positions_dynodes_bas.append((x_start, x_end, y_start, y_end))

    for i in range(N//2):
        x_start = a + (i+1)*c + d/2 + i*d - c/2
        x_end = x_start + c
        y_start = f-b
        y_end = y_start + e
        positions_dynodes_haut.append((x_start, x_end, y_start, y_end))

    for _ in range(it_max):
        Ex_val, Ey_val = Eulerer_lechamp(x[-1], y[-1], Ex, Ey, dx)

        ax = (q * Ex_val) / m
        ay = (q * Ey_val) / m

        vx += ax * dt
        vy += ay * dt

        x_old = x[-1]
        y_old = y[-1]

        x_new = x_old + vx * dt
        y_new = y_old + vy * dt

        if not (0 <= x_new < Lx and -Ly/2 <= y_new <= Ly/2):
            print("L'électron est sorti du tube")
            break

        collision_bas = False
        collision_haut = False

        for (x_start, x_end, y_start, y_end) in positions_dynodes_bas:
            if segment_intersect_rectangle(x_old, y_old, x_new, y_new, x_start, x_end, y_start, y_end):
                collision_bas = True
                break

        for (x_start, x_end, y_start, y_end) in positions_dynodes_haut:
            if segment_intersect_rectangle(x_old, y_old, x_new, y_new, x_start, x_end, y_start, y_end):
                collision_haut = True
                break

        if collision_bas:
            print(f"Collision dynode bas à ({x_new:.2f}, {y_new:.2f})")
            y_new += 2
            vy = -vy

        elif collision_haut:
            print(f"Collision dynode haut à ({x_new:.2f}, {y_new:.2f})")
            y_new -= 2
            vy = -vy

        x.append(x_new)
        y.append(y_new)

    return np.array(x), np.array(y)

# --------- SIMULATION ---------
x0 = 0.1
y0 = 0.1
vx0 = 1500
vy0 = 3
dt = 1e-5
it_max = 100000

traj_x, traj_y = position_el(x0, y0, vx0, vy0, Ex, Ey, dx, dt, it_max, a, b, c, d, e, f, N)

# --------- AFFICHAGE ---------
plt.figure(figsize=(10, 6))
cp = plt.contourf(X, Y, res, levels=100, cmap='plasma')
plt.colorbar(cp, label="Potentiel (V)")

plt.contour(X, Y, bloqué, levels=[0.5], colors='black', linewidths=0.8)

plt.plot(traj_x, traj_y, 'y-', label="Trajectoire")
plt.plot(x0, y0, 'go', label="Départ", markersize=5)

plt.xlabel("x (mm)")
plt.ylabel("y (mm)")
plt.title("Q3c - Trajectoire avec rebonds dynodes")
plt.legend()
plt.axis("equal")
plt.tight_layout()
plt.savefig("Q3c_rebonds_FINAL.png", dpi=300)
plt.show()

# --------- ANIMATION ---------

fig, ax = plt.subplots(figsize=(10, 6))
cp = ax.contourf(X, Y, res, levels=100, cmap='plasma')
plt.colorbar(cp, label="Potentiel (V)", ax=ax)
ax.contour(X, Y, bloqué, levels=[0.5], colors='black', linewidths=0.8)

electron, = ax.plot([], [], 'go', markersize=5, label="Électron")
trace, = ax.plot([], [], 'y-', linewidth=1, label="Trajectoire")

ax.set_xlabel("x (mm)")
ax.set_ylabel("y (mm)")
ax.set_title("Q3c - Animation de la trajectoire")
ax.legend()
ax.axis("equal")
ax.set_xlim(0, Lx)
ax.set_ylim(-Ly/2, Ly/2)

# Initialisation
def init():
    electron.set_data([], [])
    trace.set_data([], [])
    return electron, trace

# Animation frame
def update(frame):
    if frame < len(traj_x):
        electron.set_data(traj_x[frame], traj_y[frame])
        trace.set_data(traj_x[:frame+1], traj_y[:frame+1])
    return electron, trace

ani = animation.FuncAnimation(fig, update, frames=len(traj_x), init_func=init, blit=True, interval=1)

plt.tight_layout()
plt.show()

# Sauvegarder l'animation (optionnel, peut être lent)
# ani.save('trajectoire_electron.gif', fps=60)



print(f"Position finale : ({traj_x[-1]:.3f}, {traj_y[-1]:.3f}) mm")
print(f"Déplacement total : Δx = {traj_x[-1] - traj_x[0]:.3f} mm, Δy = {traj_y[-1] - traj_y[0]:.3f} mm")
