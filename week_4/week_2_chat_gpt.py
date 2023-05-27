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

def simulate(_f_res):
    # Array met krachten (sinusvormige oscillatie)
    F = Fmax * np.sin(2 * np.pi * _f_res * t)

    # Arrays voor verplaatsing en snelheid
    x = np.zeros_like(t)
    v = np.zeros_like(t)

    # Numerieke integratie van het Duffing-model met de Runge-Kutta-methode
    for i in range(1, len(t)):
        k1x = v[i - 1]
        k1v = (F[i - 1] - gamma * v[i - 1] - k * x[i - 1] - x[i - 1] ** 3) / m

        k2x = v[i - 1] + 0.5 * dt * k1v
        k2v = (F[i - 1] - gamma * (v[i - 1] + 0.5 * dt * k1v) - k * (x[i - 1] + 0.5 * dt * k1x) - (
                    x[i - 1] + 0.5 * dt * k1x) ** 3) / m

        k3x = v[i - 1] + 0.5 * dt * k2v
        k3v = (F[i - 1] - gamma * (v[i - 1] + 0.5 * dt * k2v) - k * (x[i - 1] + 0.5 * dt * k2x) - (
                    x[i - 1] + 0.5 * dt * k2x) ** 3) / m

        k4x = v[i - 1] + dt * k3v
        k4v = (F[i - 1] - gamma * (v[i - 1] + dt * k3v) - k * (x[i - 1] + dt * k3x) - (x[i - 1] + dt * k3x) ** 3) / m

        x[i] = x[i - 1] + (1 / 6) * dt * (k1x + 2 * k2x + 2 * k3x + k4x)
        v[i] = v[i - 1] + (1 / 6) * dt * (k1v + 2 * k2v + 2 * k3v + k4v)
    return x


# Plot de uitwijking van de massa
plt.plot(t, simulate(f_res))
plt.xlabel('Tijd (s)')
plt.ylabel('Uitwijking (m)')
plt.title('Uitwijking van de massa')
plt.show()

print('Amplitude max: ' + str(np.max(simulate(f_res))))

def simulate_frequencies(num_frequencies, frequency_step_size):
    # Frequenties rond de resonantiefrequentie
    _frequencies = np.linspace(f_res - frequency_step_size * f_res, f_res + frequency_step_size * f_res, num_frequencies)
    print(_frequencies)
    _amplitudes = np.zeros_like(_frequencies)
    for i in range(0, len(_frequencies)):
        _amplitudes[i] = np.max(simulate(_frequencies[i]))

    return _frequencies, _amplitudes


# Amplitude x, against frequencies
frequencies, amplitudes = simulate_frequencies(10000, 0.8)

# Bepaal de FWHM
half_max_amplitude = 0.5 * np.max(amplitudes)
mask = amplitudes >= half_max_amplitude
left_index = np.argmax(mask)
right_index = len(mask) - np.argmax(mask[::-1]) - 1
FWHM = frequencies[right_index] - frequencies[left_index]

# plt.plot(frequencies, amplitudes)
# plt.axhline(y=half_max_amplitude, color='r', linestyle='--', label='Half maximum')
# plt.xlabel('Frequentie (Hz)')
# plt.ylabel('Trillingsamplitude (m)')
# plt.title('Tuning-kromme van de drive mode')
# plt.legend()
# plt.text(frequencies[left_index], half_max_amplitude, f'FWHM: {FWHM:.3f} Hz', ha='left', va='bottom')
# plt.show()

# Plot the tuning curve
plt.plot(frequencies, amplitudes)
plt.xlabel('Frequentie (Hz)')
plt.ylabel('Trillingsamplitude (m)')
plt.title('Tuning-kromme van de drive mode')

# Plot horizontal line at half maximum
# plt.axhline(y=half_max_amplitude, color='r', linestyle='--', label='Half maximum')

# Find the maximum amplitude index
max_amplitude_index = np.argmax(amplitudes)
max_amplitude_frequency = frequencies[max_amplitude_index]
max_amplitude = amplitudes[max_amplitude_index]

# Plot vertical line at the maximum amplitude
# plt.axvline(x=frequencies[left_index], ymin=0, ymax=max_amplitude / np.max(amplitudes), color='g', linestyle='--', label='Max amplitude')
# plt.axvline(x=frequencies[right_index], ymin=0, ymax=max_amplitude / np.max(amplitudes), color='g', linestyle='--', label='Max amplitude')

print(frequencies[left_index])
print(frequencies[right_index])

# Plot FWHM region
# plt.fill_between(frequencies[left_index:right_index+1], amplitudes[left_index:right_index+1], alpha=0.3, color='b', label='FWHM')

# Plot legend
# plt.legend()

# Add an arrow between left and right indices
arrow_y = half_max_amplitude + 0.05 * np.max(amplitudes)
arrow_x_start = frequencies[left_index]
arrow_x_end = frequencies[right_index]
arrow_length = arrow_x_end - arrow_x_start
arrow_props = dict(arrowstyle='|-|', facecolor='black', edgecolor='black')
plt.annotate('', xy=(arrow_x_start, arrow_y), xytext=(arrow_x_end, arrow_y), arrowprops=arrow_props)

# Add text for the arrow length
arrow_text_x = arrow_x_start + 0.5 * arrow_length
plt.text(arrow_text_x + 500, arrow_y + 0.02 * np.max(amplitudes), f'FWHM: {arrow_length:.3f} Hz', ha='center', va='bottom')

plt.show()