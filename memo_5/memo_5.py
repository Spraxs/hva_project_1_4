import math

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


def simulation_timeline(m, k, Nosc, Npo):
    f = 1 / (2 * np.pi * math.sqrt(m / k))
    w = 2 * np.pi * f
    T = 1 / f
    dt = T / Npo
    t = np.arange(0, Nosc * T, dt)
    return t, dt, w, T

def simulate_voltage(t):
    Vt = np.zeros(len(t))
    for i in range(0, len(t)):
        Vt[i] = (V_0 + V_drive * np.cos(w_drive * t[i]))
    return Vt

def simulate_electrical_force(t):
    F_el = np.zeros(len(t))
    for i in range(0, len(t)):
        F_el[i] = 2 * N_drive * epsilon_zero * w * V_0 * V_drive * np.cos(w_drive * t[i]) / (d_drive)
    return F_el


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


def plot_simulation_voltage(t, voltage_array):
    # Plot de input spanning tegen tijd
    plt.plot(t, voltage_array)
    plt.legend()
    plt.xlabel('Tijd (s)')
    plt.ylabel('Voltage (V)')
    plt.title('Input Spanning')
    plt.show()


def plot_simulation_electrical_force(t, uitwijking_drive):
    # Plot de elektrische kracht
    plt.plot(t, uitwijking_drive)
    plt.legend()
    plt.xlabel('Tijd (s)')
    plt.ylabel('Elektrische kracht (N)')
    plt.title('Elektrische kracht')
    plt.show()


def plot_simulation_drive_mode(t, f_el_array):
    # Plot de uitwijking van de massa
    plt.plot(t, f_el_array)
    plt.legend()
    plt.xlabel('Tijd (s)')
    plt.ylabel('Uitwijking (m)')
    plt.title('Uitwijking (drive mode)')
    plt.show()


time_line, delta_time, w_drive, T = simulation_timeline(m_drive, k_drive, 248, 500)
voltage_array = simulate_voltage(time_line)
electrical_force_array = simulate_electrical_force(time_line)

plot_simulation_voltage(time_line, voltage_array)
plot_simulation_electrical_force(time_line, electrical_force_array)

x_drive = simulate_drive_mode_by_electrical_force(time_line, delta_time,electrical_force_array)
plot_simulation_drive_mode(time_line, x_drive)