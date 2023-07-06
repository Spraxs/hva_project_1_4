import math

import numpy as np
import matplotlib.pyplot as plt

# Stap 1 parameters:
# 0 spanning (v, volt)
V_0 = 15
# Drive spanning (v, volt)
V_drive = 1.5
# Lengte overlappend deel van condensator plaat in evenwicht (m)
l_drive = 200e-06
# Afstand tussen condensator platen (m)
d_drive = 2e-06
# Breedte van condensator platen (m) [bekend als 'h' in word notities]
w = 3e-06
# Aantal blauwe condensatorplaten
N_drive = 100
# ε0
epsilon_zero = 8.854e-12
# Elektrische veld constante
k_el = 8.988e+09

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

# Stap 5 parameters
# Lengte overlappend deel van condensator plaat in evenwicht (m)
l_sense = 200e-06
# Afstand tussen condensator platen (m)
d_sense = 2e-06
# Aantal blauwe condensatorplaten
N_sense = 40
# Spanning DC (v, volt)
Vdc = 15


def simulation_timeline(m, k, amount_of_oscillations, accuracy):
    f = 1 / (2 * np.pi * math.sqrt(m / k))
    w = 2 * np.pi * f
    T = 1 / f
    dt = T / accuracy
    t = np.arange(0, amount_of_oscillations * T, dt)
    return t, dt, w, T


def simulate_voltage(t):
    voltage_array = np.zeros(len(t))
    for i in range(0, len(t)):
        voltage_array[i] = (V_0 + V_drive * np.cos(w_drive * t[i]))
    return voltage_array


def simulate_electrical_force(t):
    F_el = np.zeros(len(t))
    for i in range(0, len(t)):
        F_el[i] = 2 * N_drive * epsilon_zero * k_el * w * V_0 * V_drive * np.cos(w_drive * t[i]) / d_drive
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


def simulate_coriolis_force(m, _x_drive):
    y_values = np.zeros(len(x_drive))
    z_values = np.zeros(len(x_drive))
    # Horizontale beweging op de x-as
    v_drive_on_x_axis = w_drive * _x_drive
    v_drive_coordinate_vector = np.column_stack((v_drive_on_x_axis, y_values, z_values))
    # Rotatie rond de z-as
    omega_gyro_vector = np.array([0, 0, omega_gyro])

    # Bereken coriolis met kruisproduct
    return np.cross(-2 * m * omega_gyro_vector, v_drive_coordinate_vector)


# Solving: 𝑚𝑥̈+ 𝛾𝑥̇ + 𝑘𝑥 = F_coriolis
def simulate_sense_mode_by_coriolis_force_y_axis(t_array, dt, f_coriolis):
    f_coriolis_y_values = f_coriolis[:, 1]

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
        force = f_coriolis_y_values[i]

        # v'(t) = ((F - c * v - k * x) / m)
        # v(t+dt) = v(t) + dt * v'(t)
        v_array[i + 1] = v + dt * ((force - gamma_sense * v - k_sense * x) / m_sense)

        # x(t+dt) = x(t) + dt * v(t)
        x_array[i + 1] = x + dt * v

    return x_array


def simulate_sense_current(t_array, x_sense):
    current_array = np.zeros(len(t_array) - 1)
    c = k_el * epsilon_zero * w / d_sense * Vdc * N_sense
    for i in range(0, len(t_array) - 1):
        q_0 = c * (x_sense[i] + l_sense)
        q_1 = c * (x_sense[i+1] + l_sense)

        t_0 = time_line[i]
        t_1 = time_line[i+1]

        dq = (q_1 - q_0)
        dt = (t_1 - t_0)

        current_array[i] = dq / dt

    return current_array


def plot_simulation_voltage(t, _voltage_array):
    plt.plot(t, _voltage_array)
    plt.xlabel('Tijd (s)')
    plt.ylabel('Voltage (V)')
    plt.title('Input Spanning tegen tijd')
    plt.show()


def plot_simulation_electrical_force(t, f_el):
    plt.plot(t, f_el)
    plt.xlabel('Tijd (s)')
    plt.ylabel('Elektrische kracht (N)')
    plt.title('Elektrische kracht tegen tijd')
    plt.show()


def plot_simulation_drive_mode(t, _x_drive):
    plt.plot(t, _x_drive)
    plt.xlabel('Tijd (s)')
    plt.ylabel('Drive uitwijking X-as (m)')
    plt.title('Drive mode uitwijking tegen tijd')
    plt.show()


def plot_coriolis_force(t, f_coriolis):
    coriolis_force_y_axis = f_coriolis[:, 1]
    plt.plot(t, coriolis_force_y_axis)
    plt.xlabel('Tijd (s)')
    plt.ylabel('F coriolis Y-axis (N)')
    plt.title('Coriolis kracht tegen tijd')
    plt.show()


def plot_x_sense(t, _x_sense_y_axis):
    plt.plot(t, _x_sense_y_axis)
    plt.xlabel('Tijd (s)')
    plt.ylabel('Sense uitwijking Y-as (m)')
    plt.title('Drive mode uitwijking tegen tijd')
    plt.show()


def plot_sense_current(t, _sense_current):
    plt.plot(t[:-1], _sense_current)
    plt.xlabel('Tijd (s)')
    plt.ylabel('Sense stroom (A)')
    plt.title('Sense stroom tegen tijd')
    plt.show()


def print_transmission_coefficients(force_electrical, input_voltage, _x_drive, force_coriolis, _x_sense, _sense_current):
    print("TC 1", np.max(force_electrical) / np.max(input_voltage))
    print("TC 2", np.max(_x_drive) / np.max(force_electrical))
    print("TC 3", np.max(force_coriolis) / np.max(_x_drive))
    print("TC 4", np.max(_x_sense) / np.max(force_coriolis))
    print("TC 5", np.max(_sense_current) / np.max(_x_sense))


time_line, delta_time, w_drive, T = simulation_timeline(m_drive, k_drive, 300, 600)

# Stap 1
voltage_array = simulate_voltage(time_line)
electrical_force_array = simulate_electrical_force(time_line)

plot_simulation_voltage(time_line, voltage_array)
plot_simulation_electrical_force(time_line, electrical_force_array)

# Stap 2
x_drive = simulate_drive_mode_by_electrical_force(time_line, delta_time,electrical_force_array)
plot_simulation_drive_mode(time_line, x_drive)

# Stap 3
coriolis_force = simulate_coriolis_force(m_drive, x_drive)
plot_coriolis_force(time_line, coriolis_force)

# Stap 4
x_sense_y_axis = simulate_sense_mode_by_coriolis_force_y_axis(time_line, delta_time, coriolis_force)
plot_x_sense(time_line, x_sense_y_axis)

# Stap 5
sense_current = simulate_sense_current(time_line, x_sense_y_axis)
plot_sense_current(time_line, sense_current)

# Print Transmissie coefficienten
print_transmission_coefficients(electrical_force_array, voltage_array, x_drive, coriolis_force, x_sense_y_axis, sense_current)