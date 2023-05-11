import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# veerconstante, k (N/m)
k = 40

# massa, m (kg)
m = 1.845 * 10 ** -6

# gamma, c (kg/s)
c = 0.042953

colNames = ['Tijd', 'Afstand']
df_xf_1 = pd.read_csv("../posities_1_Team_B6.txt", names=colNames, delimiter="   ")
df_xf_2 = pd.read_csv("../posities_2_Team_B6.txt", names=colNames, delimiter="   ")


def diff_pos_to_acc(df_xf):
    dx_array = df_xf['Afstand'].diff()
    # tijd_v array
    dt_array = df_xf['Tijd'].diff()
    # snelheid array
    velocity_array = dx_array / dt_array

    df_vf = pd.DataFrame()
    df_vf['Tijd'] = df_xf['Tijd'].iloc[1:]
    df_vf['Snelheid'] = velocity_array.iloc[1:]

    dv_array = df_vf['Snelheid'].diff()
    # tijd_a array
    dt_array = df_vf['Tijd'].diff()
    # versnelling array
    acceleration_array = dv_array / dt_array

    df_af = pd.DataFrame()
    df_af['Tijd'] = df_vf['Tijd'].iloc[1:]
    df_af['Versnelling'] = acceleration_array.iloc[1:]

    return df_af


# TODO Better function naming
# 𝑚𝑥̈+ 𝛾𝑥̇ + 𝑘𝑥 = 𝑚𝑥̈_f
def solve_for_x(df_af):
    # tijd
    time_array = df_af['Tijd'].to_numpy()
    dt = np.diff(time_array)[0]

    # begin waardes
    x0 = 0
    v0 = 0

    # positie en eerste stap
    x_array = np.zeros(len(time_array))
    v_array = np.zeros(len(time_array))

    x_array[0] = x0
    x_array[0] = v0

    for i in range(1, len(time_array) - 1):
        x = x_array[i]
        v = v_array[i]
        af = df_af['Versnelling'].iloc[i]
        # v'(t) = ((m*af - c * v - k * x) / m)
        # v(t+dt) = v(t) + dt * v'(t)
        v_array[i + 1] = v + dt * ((m * af - c * v - k * x) / m)

        # x(t+dt) = x(t) + dt * v(t)
        x_array[i + 1] = x + dt * v

    return time_array, x_array


def plot_response(time_array, x_array, name):
    # Plotting the graph
    plt.plot(time_array, x_array)

    # Adding title and labels to the graph
    plt.title(name)
    plt.xlabel("Tijd (s)")
    plt.ylabel("Respons (m)")

    # Displaying the graph
    plt.show()


def plot_response_vs_acc(time_array, x_array, af_array, name):
    fig, axs = plt.subplots(1, 2, figsize=(10, 5))
    fig.subplots_adjust(wspace=0.5)

    axs[0].plot(time_array, af_array)
    axs[0].set_title('Versnelling tegen tijd')
    axs[0].set_xlabel('Tijd (s)')
    axs[0].set_ylabel('Versnelling (m/s²)')

    axs[1].plot(time_array, x_array)
    axs[1].set_title('Respons tegen tijd')
    axs[1].set_xlabel('Tijd (s)')
    axs[1].set_ylabel('Respons (m)')

    # Adding title and labels to the graph
    plt.suptitle(name)

    # Displaying the graph
    plt.show()


df_af_1 = diff_pos_to_acc(df_xf_1)
df_af_2 = diff_pos_to_acc(df_xf_2)

time_array_1, x_array_1 = solve_for_x(df_af_1)
a_1 = df_af_1['Versnelling'].to_numpy()
scale_1 = a_1.max() / x_array_1.max()
plot_response_vs_acc(time_array_1, x_array_1 * scale_1, a_1, "Versnellingsprofiel 1")

time_array_2, x_array_2 = solve_for_x(df_af_2)
a_2 = df_af_2['Versnelling'].to_numpy()
scale_2 = a_2.max() / x_array_2.max()
plot_response_vs_acc(time_array_2, x_array_2 * scale_2, a_2, "Versnellingsprofiel 2")

