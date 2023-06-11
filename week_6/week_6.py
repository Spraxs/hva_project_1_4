import numpy as np
import matplotlib.pyplot as plt

# Frequentie van drive mode (Hz)
nu_drive = 4000
# sense mode massa (kg)
massa = 4.9e-09
# dempingscoëfficiënt van sense mode (kg/s)
gamma = 1.17e-06

# Geschaalde Coriolis kracht N s/rad
# Ook wel f_max
f_coriolis = 8e-12

k = 3.095


# Dit rotatie profiel werd gebruikt voor opdracht 4
def rotatie_snelheid_array(t):
    min_value = 2
    max_value = 5
    array = np.zeros(len(t))
    for i in range(0, len(t)):
        if i < len(t) / 8 or i > len(t) * 5 / 8:
            array[i] = min_value
        else:
            array[i] = max_value
    return array


# Dit rotatie profiel werd gebruikt voor de respons tijd te berekenen voor opdracht 8
def rotatie_snelheid_respons_tijd_array(t):
    min_value = 0
    max_value = 0.1
    array = np.zeros(len(t))
    for i in range(0, len(t)):
        if i < len(t) / 8 * 2:
            array[i] = min_value
        else:
            array[i] = max_value
    return array


# resonantie frequentie = drive frequentie
def force_coriolis_function(t, _rotatie_snelheid_arr):
    # (rad/s)
    hoek_snelheid = 2 * np.pi * nu_drive
    array = np.zeros(len(t))
    for i in range(0, len(t)):
        # Ω (rad / s)
        array[i] = f_coriolis * _rotatie_snelheid_arr[i] * np.cos(hoek_snelheid * t[i])
    return array


def simulation_timeline():
    # Totale simulatieduur (s)
    t_max = 0.2
    # Tijdstap (aantal tijdstappen aanpassen voor gewenste nauwkeurigheid)
    delta_time = 0.000000007
    return np.arange(0, t_max, delta_time), delta_time


# Solving: 𝑚𝑥̈+ 𝛾𝑥̇ + 𝑘𝑥 = F_coriolis
def solve_differential_equation(t_array, delta_time, force_coriolis_array, _gamma):
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
        v_array[i + 1] = v + delta_time * ((force - _gamma * v - k * x) / massa)

        # x(t+dt) = x(t) + dt * v(t)
        x_array[i + 1] = x + delta_time * v

    return x_array


def plot_simulation(t, normal_gamma_values, small_gamma_values, big_gamma_values):
    print("Totaal tijd: " + str(t[-1]))
    # Plot de uitwijking van de massa
    plt.plot(t, small_gamma_values, label="Kleine gamma")
    plt.plot(t, normal_gamma_values, label="Normale gamma")
    plt.plot(t, big_gamma_values, label="Grote gamma")
    plt.legend()
    plt.xlabel('Tijd (s)')
    plt.ylabel('Uitwijking (m)')
    plt.title('Uitwijking van de massa')
    plt.show()


def plot_f(t, kracht):
    print("Totaal tijd: " + str(t[-1]))
    # Plot de uitwijking van de massa
    plt.plot(t, kracht)
    plt.xlabel('Tijd (s)')
    plt.ylabel('Kracht (m)')
    plt.title('Uitwijking van de massa')
    plt.show()


def plot_rotatiesnelheid(t, rotatie_snelheid):
    plt.plot(t, rotatie_snelheid)
    plt.xlabel('Tijd (s)')
    plt.ylabel('Rotatiesnelheid (rad/s)')
    plt.title('Rotatie profiel')
    plt.show()


time_array, dt = simulation_timeline()
rotatie_snelheid_arr = rotatie_snelheid_respons_tijd_array(time_array)
plot_rotatiesnelheid(time_array, rotatie_snelheid_arr)
coriolis_force_array = force_coriolis_function(time_array, rotatie_snelheid_arr)

solved_differential_normal_gamma_array = solve_differential_equation(time_array, dt, coriolis_force_array, gamma)
solved_differential_small_gamma_array = solve_differential_equation(time_array, dt, coriolis_force_array, gamma / 2)
solved_differential_big_gamma_array = solve_differential_equation(time_array, dt, coriolis_force_array, gamma * 2)

plot_simulation(
    t=time_array,
    normal_gamma_values=solved_differential_normal_gamma_array,
    small_gamma_values=solved_differential_small_gamma_array,
    big_gamma_values=solved_differential_big_gamma_array
)
