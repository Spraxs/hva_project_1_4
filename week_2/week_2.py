import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import result_holder as result

colNames = ['Tijd', 'Afstand']
df_xf_1 = pd.read_csv("../posities_1_Team_B6.txt", names=colNames, delimiter="   ")
df_xf_2 = pd.read_csv("../posities_2_Team_B6.txt", names=colNames, delimiter="   ")

# veerconstante, k (N/m)
# massa, m (kg)
# gamma, c (kg/s)
k = 40
k_edited = 45
m = 1.845 * 10 ** -6
c = 0.042953
kritisch_gedempt_c = 1.59 * 10 ** (-8)


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


# Solving: 𝑚𝑥̈+ 𝛾𝑥̇ + 𝑘𝑥 = 𝑚𝑥̈_f
def simulate_accelerometer(time_array, af_array, _c, _m, _k):
    dt = np.diff(time_array)[0]

    # Boundary conditions
    x0 = 0
    v0 = 0

    # Array with boundary conditions
    x_array = np.zeros(len(time_array))
    v_array = np.zeros(len(time_array))

    x_array[0] = x0
    x_array[0] = v0

    # Solve numerically
    for i in range(0, len(time_array) - 1):
        x = x_array[i]
        v = v_array[i]
        af = af_array[i]
        # v'(t) = ((m*af - c * v - k * x) / m)
        # v(t+dt) = v(t) + dt * v'(t)
        v_array[i + 1] = v + dt * ((- _c * v - _k * x + _m * af) / _m)

        # x(t+dt) = x(t) + dt * v(t)
        x_array[i + 1] = x + dt * v

    return x_array


def plot_response(time_array, x_array, name):
    # Plotting the graph
    plt.plot(time_array, x_array)

    # Adding title and labels to the graph
    plt.title(name)
    plt.xlabel("Tijd (s)")
    plt.ylabel("Respons (m)")

    # Displaying the graph
    plt.show()


def plot_acc_with_ideal_response(_results, name):
    fig, axs = plt.subplots(len(results), 1, figsize=(10, 5))
    fig.subplots_adjust(hspace=1)
    for i in range(0, len(_results)):
        r = _results[i]
        response_ideal = m * r.af_array / k

        axs[i].plot(r.time_array, r.response_array, label='Originele respons')
        # axs[i].plot(r.time_array, r.response_edited, label='Kritisch gedempte respons')
        axs[i].plot(r.time_array, r.response_edited_with_k, label='k & gamma respons')
        axs[i].plot(r.time_array, response_ideal, label='Ideale respons')

        axs[i].set_title('Profiel ' + str(i + 1))
        axs[i].set_xlabel('Tijd (s)')
        axs[i].set_ylabel('x (m)')
        axs[i].legend()

    # Adding title and labels to the graph
    plt.suptitle(name)

    # Displaying the graph
    plt.show()


def calculate_result(df_xf):
    df_vf, df_af = diff_pos_to_acc(df_xf)

    time_array = df_af['Tijd'].to_numpy()
    af_array = df_af['Versnelling'].to_numpy()

    response_array = simulate_accelerometer(time_array, af_array, c, m, k)
    response_array_edited = simulate_accelerometer(time_array, af_array, kritisch_gedempt_c, m, k)
    response_array_edited_with_k = simulate_accelerometer(time_array, af_array, kritisch_gedempt_c, m, k_edited)

    return result.Holder(time_array, response_array, response_array_edited, response_array_edited_with_k, af_array)


result_1 = calculate_result(df_xf_1)
result_2 = calculate_result(df_xf_2)
results = [result_1, result_2]

plot_acc_with_ideal_response(results, 'Responsen tegen tijd')
