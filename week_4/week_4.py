import numpy as np


def force(t_array, frequency):
    return np.sin(2*np.pi*frequency*t_array)


# veerconstante, k (N/m)
# massa, m (kg)
# gamma, c (kg/s)
# amplitude, F_max (N)
k = 0.2013
m = 5.1 * 10 ** -9
c = 4.9299 * 10 ** -7
F_max = 60 * 10 ** -9

f = 1 / (2 * np.pi * np.sqrt(m/k))

t_max = 100
N_steps = 25+1
time_array = np.linspace(0, t_max, N_steps)
print(time_array)


F_array = force(time_array, f)
print(F_array)


