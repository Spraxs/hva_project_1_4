import math

import numpy as np
import matplotlib.pyplot as plt

k = 40 # N/m
m = 0.001845 # kg
gamma = 0.042953 # kg/s

A_0 = 1

t_start = 0
t_end = 10
dt = 0.01
t = np.arange(t_start, t_end, dt)

def offset_x(t):
    b = 2 * gamma * m
    result = math.sqrt(4 * m * k)
    w = math.sqrt((k/m) - b**2/(4*m**2))
    return A_0 * math.e**(-gamma*t) * math.cos(w*t)

def position(t):
    x = np.zeros_like(t)
    v = np.zeros_like(t)
    for i in range(1, len(t)):
        v[i] = v[i-1] + (dt/m) * (-k*x[i-1] - gamma*v[i-1] + acceleration(t[i-1]))
        x[i] = x[i-1] + dt * v[i]
    return x

fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, figsize=(10, 8))
fig.suptitle('Accelerometer simulatie')
ax1.plot(t, acceleration(t), 'g')
ax1.set_ylabel('Versnelling (m/s^2)')
ax1.set_xlabel('Tijd (s)')
ax2.plot(t, position(t), 'b')
ax2.set_ylabel('Positie (m)')
ax2.set_xlabel('Tijd (s)')
ax3.plot(t, velocity(t), 'r')
ax3.set_ylabel('Snelheid (m/s)')
ax3.set_xlabel('Tijd (s)')
plt.show()
