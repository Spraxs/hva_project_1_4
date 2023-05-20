import numpy as np
import matplotlib.pyplot as plt


# Solving: 𝑚𝑥̈+ 𝛾𝑥̇ + 𝑘𝑥 = F_max * sin(2πft)
def simulate_accelerometer(t_array, frequency, force_max, _c, _m, _k, delta_time):
    # Boundary conditions
    x0 = 0
    v0 = 0

    # Array with boundary conditions
    x_array = np.zeros(len(t_array))
    v_array = np.zeros(len(t_array))

    x_array[0] = x0
    x_array[0] = v0

    # Solve numerically
    for i in range(0, len(t_array) - 1):
        x = x_array[i]
        v = v_array[i]
        t = t_array[i]
        # v'(t) = ((F - c * v - k * x) / m)
        # v(t+dt) = v(t) + dt * v'(t)
        v_array[i + 1] = v + delta_time * ((force_max * np.sin(2*np.pi*frequency*t) - _c * v - _k * x) / _m)

        # x(t+dt) = x(t) + dt * v(t)
        x_array[i + 1] = x + delta_time * v

    return x_array


# veerconstante, k (N/m)
# massa, m (kg)
# gamma, c (kg/s)
# amplitude, F_max (N)
k = 0.2013
m = 5.1 * 10 ** -9
c = 4.9299 * 10 ** -7
F_max = 60 * 10 ** -9

f = 1 / (2 * np.pi * np.sqrt(m/k))

t_max = 0.005
N_steps = 100+1
time_array = np.linspace(0, t_max, N_steps)
dt = t_max/(N_steps-1)

print(time_array)

response_array = simulate_accelerometer(time_array, f, F_max, c, m, k, dt)

# Plotting the graph
plt.plot(time_array, response_array)

# Adding labels and title
plt.xlabel('Tijd (s)')
plt.ylabel('Respons (m)')
plt.title('Simulatie')

# Display the graph
plt.show()

