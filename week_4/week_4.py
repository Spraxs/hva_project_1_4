import numpy as np
import matplotlib.pyplot as plt

# Parameters
m = 5.1e-9  # massa in kg
k = 0.2013  # veerconstante in N/m
gamma = 4.9299e-7  # dempingsconstante in kg/s
Fmax = 60e-9  # amplitude van de kracht in N
f_res = 1 / (2 * np.pi * np.sqrt(m / k))  # resonantiefrequentie

# Simulatie instellingen
tmax = 90 * (1 / f_res)  # Totale simulatieduur
dt = tmax / 500  # Tijdstap (aantal tijdstappen aanpassen voor gewenste nauwkeurigheid)

# Array met tijdstappen
t = np.arange(0, tmax, dt)

# gamma += gamma * 0.2

def simulate(_f_res, _gamma):
    # Array met krachten (sinusvormige oscillatie)
    F = Fmax * np.sin(2 * np.pi * _f_res * t)

    # Arrays voor verplaatsing en snelheid
    x = np.zeros_like(t)
    v = np.zeros_like(t)

    # Numerieke integratie van het Duffing-model met de Runge-Kutta-methode
    for i in range(1, len(t)):
        k1x = v[i - 1]
        k1v = (F[i - 1] - _gamma * v[i - 1] - k * x[i - 1] - x[i - 1] ** 3) / m

        k2x = v[i - 1] + 0.5 * dt * k1v
        k2v = (F[i - 1] - _gamma * (v[i - 1] + 0.5 * dt * k1v) - k * (x[i - 1] + 0.5 * dt * k1x) - (
                    x[i - 1] + 0.5 * dt * k1x) ** 3) / m

        k3x = v[i - 1] + 0.5 * dt * k2v
        k3v = (F[i - 1] - _gamma * (v[i - 1] + 0.5 * dt * k2v) - k * (x[i - 1] + 0.5 * dt * k2x) - (
                    x[i - 1] + 0.5 * dt * k2x) ** 3) / m

        k4x = v[i - 1] + dt * k3v
        k4v = (F[i - 1] - _gamma * (v[i - 1] + dt * k3v) - k * (x[i - 1] + dt * k3x) - (x[i - 1] + dt * k3x) ** 3) / m

        x[i] = x[i - 1] + (1 / 6) * dt * (k1x + 2 * k2x + 2 * k3x + k4x)
        v[i] = v[i - 1] + (1 / 6) * dt * (k1v + 2 * k2v + 2 * k3v + k4v)
    return x


def plot_simulation():
    # Plot de uitwijking van de massa
    plt.plot(t, simulate(f_res, gamma))
    plt.xlabel('Tijd (s)')
    plt.ylabel('Uitwijking (m)')
    plt.title('Uitwijking van de massa')
    plt.show()


def simulate_frequencies(_frequencies, _gamma):
    _amplitudes = np.zeros_like(_frequencies)
    for i in range(0, len(_frequencies)):
        _amplitudes[i] = np.max(simulate(_frequencies[i], _gamma))
    return _amplitudes


def calculate_FWHM(_amplitudes, _frequencies):
    # Bepaal de FWHM
    half_max_amplitude = 0.5 * np.max(_amplitudes)
    mask = _amplitudes >= half_max_amplitude
    left_index = np.argmax(mask)
    right_index = len(mask) - np.argmax(mask[::-1]) - 1
    return _frequencies[right_index] - _frequencies[left_index]


def plot_amplitudes_against_frequencies(num_frequencies, frequency_step_size, _gamma, percentage_offset):
    _frequencies = np.linspace(f_res - frequency_step_size * f_res, f_res + frequency_step_size * f_res, num_frequencies)

    gamma_offset = _gamma * percentage_offset

    amplitudes = simulate_frequencies(_frequencies, _gamma)
    amplitudes_min = simulate_frequencies(_frequencies, _gamma - gamma_offset)
    amplitudes_max = simulate_frequencies(_frequencies, _gamma + gamma_offset)

    # Plot the tuning curve
    plt.plot(_frequencies, amplitudes_min, label='-20% gamma', color='blue', linewidth=1)
    plt.plot(_frequencies, amplitudes, label='Origineel', color='orange', linewidth=1)
    plt.plot(_frequencies, amplitudes_max, label='+20% gamma', color='green', linewidth=1)
    plt.xlabel('Frequentie (Hz)')
    plt.ylabel('Trillingsamplitude (m)')
    plt.title('Tuning-kromme van de drive mode')

    max_amplitude_min = np.max(amplitudes_min)
    max_amplitude = np.max(amplitudes)
    max_amplitude_max = np.max(amplitudes_max)

    formatted_amplitude_min = '{:.2e}'.format(max_amplitude_min)
    formatted_amplitude = '{:.2e}'.format(max_amplitude)
    formatted_amplitude_max = '{:.2e}'.format(max_amplitude_max)

    plt.axhline(y=max_amplitude_min, linestyle='--', color='blue')
    plt.axhline(y=max_amplitude, linestyle='--', color='orange')
    plt.axhline(y=max_amplitude_max, linestyle='--', color='green')

    FWHM_min = calculate_FWHM(amplitudes_min, _frequencies)
    FWHM = calculate_FWHM(amplitudes, _frequencies)
    FWHM_max = calculate_FWHM(amplitudes_max, _frequencies)

    print(f'FWHM (-20% gamma): {FWHM_min:.3f} Hz')
    print(f'FWHM (origineel): {FWHM:.3f} Hz')
    print(f'FWHM (+20% gamma): {FWHM_max:.3f} Hz')

    print(f'Max amplitude (-20% gamma): {formatted_amplitude_min} m')
    print(f'Max amplitude (origineel): {formatted_amplitude} m')
    print(f'Max amplitude (+20% gamma): {formatted_amplitude_max} m')

    # Plot legend
    plt.legend(fontsize='small')

    plt.show()


plot_simulation()
plot_amplitudes_against_frequencies(num_frequencies=10000, frequency_step_size=0.8, _gamma=gamma, percentage_offset=0.2)