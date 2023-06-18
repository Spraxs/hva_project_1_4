import numpy as np
import matplotlib.pyplot as plt

# Stap 1 parameters:
# 0 spanning (v, volt)
V_0 = 15
# Drive spanning (v, volt)
V_drive = 1.5
# Lengte overlappend deel van condensator plaat in evenwicht (
l_drive = 200e-06
# Afstand tussen condensator platen (m)
d_drive = 2e-06
# Breedte van condensator platen (m) [bekend als 'h' in word notities]
w = 3e-06
# N coefficient (omschrijving van deze coefficient is aangegeven in de opdrachten)
N_drive = 100
# ε0
epsilon_zero = 8.854e-12

# Frequentie (Hz)
v_drive = 31.6

# Stap 2 parameters:
# Massa drive (kg)
m_drive = 5.1e-09
# Veerconstante drive (N/m)
k_drive = 0.2013
# Dempingsconstante drive (kg/s)
gamma_drive = 4.9299e-07

# Stap 3 parameters:
# Rotatie snelheid (rad / s)
omega_gyro = 10 / 60

# Stap 4 parameters:
# Massa sense (kg)
m_sense = 6.12e-09
# Veerconstante sense (N/m)
k_sense = 0.9664
# Dempingscoefficient sense (kg/s)
gamma_sense = 1.538e-06


def simulation_timeline():
    # Totale simulatieduur (s)
    t_max = 0.2
    # Tijdstap (aantal tijdstappen aanpassen voor gewenste nauwkeurigheid)
    dt = 0.000000007
    return np.arange(0, t_max, dt), dt


def simulate_electrical_force(t):
    array = np.zeros(len(t))
    x_constant = epsilon_zero * k_drive * w / (2*d_drive)
    drive_rotation_speed = 2 * np.pi * v_drive
    for i in range(0, len(t)):
        # Ω (rad / s)
        array[i] = x_constant * (3*V_0**2/2 + 2*V_0*V_drive*np.cos(drive_rotation_speed * t[i]) +
                                 V_0**2/2 * np.cos(2*drive_rotation_speed * t[i]))
    return array


# Solving: 𝑚𝑥̈+ 𝛾𝑥̇ + 𝑘𝑥 = F_el
def simulate_drive_mode_by_electrical_force(t_array, dt, force_el_array):
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
        force = force_el_array[i]

        # v'(t) = ((F - c * v - k * x) / m)
        # v(t+dt) = v(t) + dt * v'(t)
        v_array[i + 1] = v + dt * ((force - gamma_drive * v - k_drive * x) / m_drive)

        # x(t+dt) = x(t) + dt * v(t)
        x_array[i + 1] = x + dt * v

    return x_array


def calculate_coriolis_force(m, _x_drive):
    # Differentieer x_drive om v_drive te verkrijgen
    _v_drive = np.diff(_x_drive) / delta_time
    # Voeg een extra punt toe voor de laatste snelheid (om lengte te behouden)
    _v_drive = np.append(_v_drive, _v_drive[-1])

    # Bereken de Coriolis-kracht: F_coriolis = -2 * m * omega x v
    return -2 * m * omega_gyro * _v_drive


# Solving: 𝑚𝑥̈+ 𝛾𝑥̇ + 𝑘𝑥 = F_el
def simulate_sense_mode_by_coriolis_force(t_array, dt, force_coriolis_array):
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
        v_array[i + 1] = v + dt * ((force - gamma_sense * v - k_sense * x) / m_sense)

        # x(t+dt) = x(t) + dt * v(t)
        x_array[i + 1] = x + dt * v

    return x_array


def plot_simulation_uitwijking_drive(t, uitwijking):
    # Plot de uitwijking van de massa
    plt.plot(t, uitwijking)
    plt.legend()
    plt.xlabel('Tijd (s)')
    plt.ylabel('Uitwijking (m)')
    plt.title('Uitwijking van de massa')
    plt.show()


def plot_simulation_coriolis_force(t, f_coriolis):
    # Plot de uitwijking van de massa
    plt.plot(t, f_coriolis)
    plt.legend()
    plt.xlabel('Tijd (s)')
    plt.ylabel('Coriolis kracht (N)')
    plt.title('Coriolis kracht')
    plt.show()


def plot_simulation_uitwijking_sense(t, uitwijking_sense):
    # Plot de uitwijking van de massa
    plt.plot(t, uitwijking_sense)
    plt.legend()
    plt.xlabel('Tijd (s)')
    plt.ylabel('Uitwijking (m)')
    plt.title('Uitwijking sense')
    plt.show()


time_line, delta_time = simulation_timeline()
electrical_force_array = simulate_electrical_force(time_line)
x_drive = simulate_drive_mode_by_electrical_force(time_line, delta_time, electrical_force_array)
coriolis_force = calculate_coriolis_force(m_drive, x_drive)
x_sense = simulate_sense_mode_by_coriolis_force(time_line, delta_time, coriolis_force)

plot_simulation_uitwijking_drive(time_line, x_drive)
plot_simulation_coriolis_force(time_line, coriolis_force)
plot_simulation_uitwijking_sense(time_line, x_sense)
