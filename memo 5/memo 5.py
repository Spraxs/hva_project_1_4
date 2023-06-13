import numpy as np

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

# Stap 2 parameters:
# Massa (kg)
m_drive = 5.1e-09
# Veerconstante (N/m)
k_drive = 0.2013
# Dempingsconstante (kg/s)
gamma_drive = 4.9299e-07


def simulation_timeline():
    # Totale simulatieduur (s)
    t_max = 0.2
    # Tijdstap (aantal tijdstappen aanpassen voor gewenste nauwkeurigheid)
    delta_time = 0.000000007
    return np.arange(0, t_max, delta_time), delta_time


def simulate_electrical_force(t):
    array = np.zeros(len(t))
    x_constant = epsilon_zero * k_drive * w / (2*d_drive)
    drive_rotation_speed = 2 * np.pi *
    for i in range(0, len(t)):
        # Ω (rad / s)
        array[i] = x_constant * (3*V_0^2/2 + 2*V_0*V_drive*np.cos())
    return array

# Solving: 𝑚𝑥̈+ 𝛾𝑥̇ + 𝑘𝑥 = F_el
def simulate_drive_mode_by_electrical_force(t_array, delta_time, force_el_array):
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
        force = force_el_array[i]
        # v'(t) = ((F - c * v - k * x) / m)
        # v(t+dt) = v(t) + dt * v'(t)
        v_array[i + 1] = v + delta_time * ((force - gamma_drive * v - k_drive * x) / m_drive)

        # x(t+dt) = x(t) + dt * v(t)
        x_array[i + 1] = x + delta_time * v

    return x_array

