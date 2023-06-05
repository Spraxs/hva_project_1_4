import numpy as np
import matplotlib.pyplot as plt

# Frequentie van drive mode (Hz)
nu_drive = 4000
# sense mode massa (mg)
massa = 4.9e-06
# dempingscoëfficiënt van sense mode (mg/s)
gamma = 1.17
# Geschaalde Coriolis kracht mN s/rad
f_coriolis = 8e-09

k = 3.095e+03


# resonantie frequentie = drive frequentie

def force_coriolis_function(t):
    # Ω (rad / s)
    rotatie_snelheid = 10 / 60
    # (rad/s)
    hoek_snelheid = 1
    return f_coriolis * rotatie_snelheid * np.cos(hoek_snelheid * t)


def simulation_timeline():
    # Totale simulatieduur (s)
    t_max = 10
    # Tijdstap (aantal tijdstappen aanpassen voor gewenste nauwkeurigheid)
    delta_time = 0.01
    return np.arange(0, t_max, delta_time), delta_time


# Solving: 𝑚𝑥̈+ 𝛾𝑥̇ + 𝑘𝑥 = F_coriolis
def solve_differential_equation(t_array, delta_time, force_coriolis_array):
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
        force = force_coriolis_array[i]
        # v'(t) = ((F - c * v - k * x) / m)
        # v(t+dt) = v(t) + dt * v'(t)
        v_array[i + 1] = v + delta_time * ((force - gamma * v - k * x) / massa)

        # x(t+dt) = x(t) + dt * v(t)
        x_array[i + 1] = x + delta_time * v

    return x_array


def plot_simulation(t, uitwijking):
    print("Totaal tijd: " + str(t[-1]))
    # Plot de uitwijking van de massa
    plt.plot(t, uitwijking)
    plt.xlabel('Tijd (s)')
    plt.ylabel('Uitwijking (m)')
    plt.xlim(2.5, 5)
    plt.title('Uitwijking van de massa')
    plt.show()


time_array, dt = simulation_timeline()
coriolis_force_array = force_coriolis_function(time_array)

solved_differential_array = solve_differential_equation(time_array, dt, coriolis_force_array)


print("Time: " + str(len(time_array)) + " vs " + str(len(solved_differential_array)))
# plot_simulation(t=time_array, uitwijking=solved_differential_array)

print(solved_differential_array)


# Get the time array and delta time
time_array, delta_time = simulation_timeline()

# Example data
data = np.sin(2 * np.pi * 10 * time_array)

# Plotting
fig, ax = plt.subplots()
ax.plot(time_array, solved_differential_array)

plt.show()