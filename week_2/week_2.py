import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import ResultHolder as result

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

    return df_vf, df_af


# TODO Better function naming
# 𝑚𝑥̈+ 𝛾𝑥̇ + 𝑘𝑥 = 𝑚𝑥̈_f
def solve_for_x(df_af, c, m, k):
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
        v_array[i + 1] = v + dt * ((- c * v - k * x + m * af) / m)

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


def plot_acc_with_ideal_response(time_array, _results, name):
    fig, axs = plt.subplots(len(results), 1, figsize=(10, 5))
    fig.subplots_adjust(hspace=1)
    for i in range(0, len(_results)):
        r = _results[i]
        axs[i].plot(time_array, r.response_array, label='Originele respons')
        # TODO Un-comment:
        # axs[i].plot(time_array, r.response_edited, label='Kritisch gedempte respons')
        axs[i].plot(time_array, r.response_edited_with_k, label='k & gamma respons')
        axs[i].plot(time_array, m * r.af_array / k, label='Ideale respons')
        axs[i].set_title('Profiel ' + str(i + 1))
        axs[i].set_xlabel('Tijd (s)')
        axs[i].set_ylabel('x (m)')
        axs[i].legend()

    # Adding title and labels to the graph
    plt.suptitle(name)

    # Displaying the graph
    plt.show()

# veerconstante, k (N/m)
# massa, m (kg)
# gamma, c (kg/s)
c=0.042953
kritisch_gedempt_c = 1.59*10**(-8)
m=1.845 * 10 ** -6
k=40
k_edited=45

# TODO Clean redundant variable (like time_array_1)
df_vf_1, df_af_1 = diff_pos_to_acc(df_xf_1)
time_array_1, response_array = solve_for_x(df_af_1, c=c, m=m, k=k)
time_array_1, response_array_edited = solve_for_x(df_af_1, c=kritisch_gedempt_c, m=m, k=k)
time_array_1, response_array_edited_with_k = solve_for_x(df_af_1, c=kritisch_gedempt_c, m=m, k=k_edited)
af_array = df_af_1['Versnelling'].to_numpy()

result_1 = result.Holder(response_array, response_array_edited, response_array_edited_with_k, af_array)

df_vf_2, df_af_2 = diff_pos_to_acc(df_xf_2)
time_array_1, response_array = solve_for_x(df_af_2, c=c, m=m, k=k)
time_array_1, response_array_edited = solve_for_x(df_af_2, c=kritisch_gedempt_c, m=m, k=k)
time_array_1, response_array_edited_with_k = solve_for_x(df_af_2, c=kritisch_gedempt_c, m=m, k=k_edited)
af_array = df_af_2['Versnelling'].to_numpy()
result_2 = result.Holder(response_array, response_array_edited, response_array_edited_with_k, af_array)

results = [result_1, result_2]

plot_acc_with_ideal_response(time_array_1, results, 'Respons & ideale respons tegen tijd')