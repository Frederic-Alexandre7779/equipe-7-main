import numpy as np
import question3c

dx = 0.1
c = question3c.relaxation(V, bloqué, variation=1e-3, max_iter=3000)

# calculer le gradient en x et y grâce à numpy
Ey, Ex = np.gradient(-c, dx, dx)

E_norm = np.sqrt(Ex**2 + Ey**2) #norme

x = [2]
y = [0]

Ey_val, Ex_val = question3c.Eulerer_lechamp(x[0], y[0], Ex, Ey, dx)

print(Ex_val, Ey_val)