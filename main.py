import pandas as pd
import matplotlib.pyplot as plt

colNames = ['Tijd', 'Afstand']
df_x_1 = pd.read_csv("posities_1_Team_B6.txt", names=colNames, delimiter="   ")
df_x_2 = pd.read_csv("posities_2_Team_B6.txt", names=colNames, delimiter="   ")


def plot(df_x, name):
    fig, axs = plt.subplots(1, 3, figsize=(10, 5))
    fig.subplots_adjust(wspace=0.5)
    axs[0].plot(df_x['Tijd'], df_x['Afstand'])
    axs[0].set_title('Afstand tegen tijd ')
    axs[0].set_xlabel('Tijd')
    axs[0].set_ylabel('Afstand')

    dx_array = df_x['Afstand'].diff()
    # tijd_v array
    dt_array = df_x['Tijd'].diff()
    # snelheid array
    velocity_array = dx_array / dt_array

    df_v = pd.DataFrame()
    df_v['Tijd'] = df_x['Tijd'].iloc[1:]
    df_v['Snelheid'] = velocity_array.iloc[1:]

    axs[1].plot(df_v['Tijd'], df_v['Snelheid'])
    axs[1].set_title('Snelheid tegen tijd ')
    axs[1].set_xlabel('Tijd')
    axs[1].set_ylabel('Snelheid')

    dv_array = df_v['Snelheid'].diff()
    # tijd_a array
    dt_array = df_v['Tijd'].diff()
    # versnelling array
    acceleration_array = dv_array / dt_array

    df_a = pd.DataFrame()
    df_a['Tijd'] = df_v['Tijd'].iloc[1:]
    df_a['Versnelling'] = acceleration_array.iloc[1:]

    axs[2].plot(df_a['Tijd'], df_a['Versnelling'])
    axs[2].set_title('Versnelling tegen tijd ')
    axs[2].set_xlabel('Tijd')
    axs[2].set_ylabel('Versnelling')

    plt.suptitle(name)

    plt.show()

    return df_a['Versnelling']

# 1 array met tijd
time_values = df_x_1['Tijd'].iloc[2:]

# 2 arrays met versnelling
acceleration_values_1 = plot(df_x_1, 'Bestand 1')
acceleration_values_2 = plot(df_x_2, 'Bestand 2')

df = pd.DataFrame({'Tijd': time_values, 'Versnelling 1': acceleration_values_1, 'Versnelling 2': acceleration_values_2})

df.to_csv('versnellingen_Team_B6.csv', index=False)
