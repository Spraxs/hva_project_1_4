import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

colNames = ['Tijd', 'Afstand']
df_xf_1 = pd.read_csv("../posities_1_Team_B6.txt", names=colNames, delimiter="   ")
df_xf_2 = pd.read_csv("../posities_2_Team_B6.txt", names=colNames, delimiter="   ")

# veerconstante, k (N/m)
k = 40

m = 1.845 * 10**-6

gamma_kritisch = 1.59*10**-8

gamma_verhoogd = 1.4 * 10**2



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

    return df_vf, df_af

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

    for i in range(0, len(time_array) - 1):
        x = x_array[i]
        v = v_array[i]
        af = df_af['Versnelling'].iloc[i]
        # v'(t) = ((m*af - c * v - k * x) / m)
        # v(t+dt) = v(t) + dt * v'(t)
        v_array[i + 1] = v + dt * ((- gamma_kritisch * v - k * x + m * af) / m)

        # x(t+dt) = x(t) + dt * v(t)
        x_array[i + 1] = x + dt * v

    return time_array, x_array

def plot_acc_with_scaled_response(time_array, response_array_1, name):
    plt.plot(time_array, response_array_1)
    plt.xlabel('Tijd (s)')
    plt.ylabel('Respons (x)')

    # Adding title and labels to the graph
    plt.title(name)

    # Displaying the graph
    plt.show()

df_vf_1, df_af_1 = diff_pos_to_acc(df_xf_1)

time_array_1, response_array_1 = solve_for_x(df_af_1)
a_1 = df_af_1['Versnelling'].to_numpy()
scale_1 = a_1.max() / response_array_1.max()

plot_acc_with_scaled_response(time_array_1, response_array_1, 'Respons')