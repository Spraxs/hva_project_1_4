import numpy as np
import matplotlib.pyplot as plt


# Solving: 𝑚𝑥̈+ 𝛾𝑥̇ + 𝑘𝑥 = F_max
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
        force = force_max * np.sin(2*np.pi*frequency*t)
        v_array[i + 1] = v + delta_time * ((force - _c * v - _k * x - x**3) / _m)

        # x(t+dt) = x(t) + dt * v(t)
        x_array[i + 1] = x + delta_time * v

    return x_array


# veerconstante, k (N/m)
# massa, m (kg)
# gamma, c (kg/s)
# amplitude, F_max (N)
# m = 5.1e-11  # massa in kg
# k = 0.2013  # veerconstante in N/m
# gamma = 4.9299e-7  # dempingsconstante in kg/s
# F_max = 60e-9  # amplitude van de kracht in N

k = 0.2013
m = 5.1 * 10 ** -9
gamma = 4.9299 * 10 ** -7
F_max = 60 * 10 ** -9

f_res = 1 / (2 * np.pi * np.sqrt(m / k))  # resonantiefrequentie

# Simulatie instellingen
tmax = 50 * (1 / f_res)  # Totale simulatieduur
dt = tmax / 3000  # Tijdstap (aantal tijdstappen aanpassen voor gewenste nauwkeurigheid)

# Array met tijdstappen
time_array = np.arange(0, tmax, dt)

print(time_array)

response_array = simulate_accelerometer(time_array, f_res, F_max, gamma, m, k, dt)

# Plotting the graph
plt.plot(time_array, response_array)

# Adding labels and title
plt.xlabel('Tijd (s)')
plt.ylabel('Respons (m)')
plt.title('Simulatie')

# Display the graph
plt.show()

