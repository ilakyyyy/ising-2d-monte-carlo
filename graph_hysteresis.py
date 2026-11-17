import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

L = 128
n_runs = 15
all_data = []

for i in range(n_runs):
    filename = f"results/hysteresis_L{L}_run_{i}.txt"
    df = pd.read_csv(filename, sep='\t', names=['beta', 'en', 'mag', 'en2', 'mag2'])
    all_data.append(df)

# Split each sweep into cooling (forward) and heating (return) branches
half = len(all_data[0]) // 2
ida_list = [d.iloc[:half] for d in all_data]
vuelta_list = [d.iloc[half:] for d in all_data]

def get_stats(list_of_dfs):
    combined = pd.concat(list_of_dfs)
    mean = combined.groupby('beta').mean().reset_index()
    std = combined.groupby('beta').std().reset_index()
    return mean, std

ida_mean, ida_std = get_stats(ida_list)
vuelta_mean, vuelta_std = get_stats(vuelta_list)

plt.figure(figsize=(10, 6))
plt.errorbar(ida_mean['beta'], ida_mean['mag'], yerr=ida_std['mag'] / np.sqrt(n_runs),
             label='Enfriamiento (Ida)', fmt='-o', capsize=3)
plt.errorbar(vuelta_mean['beta'], vuelta_mean['mag'], yerr=vuelta_std['mag'] / np.sqrt(n_runs),
             label='Calentamiento (Vuelta)', fmt='-o', capsize=3)

plt.title(f"Ciclo de Histéresis Modelo de Ising (L={L}, n={n_runs})")
plt.xlabel(r"$\beta$ (1/T)")
plt.ylabel("Magnetización <|M|>")
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig(f"plots/hysteresis_L{L}.png", dpi=300, bbox_inches='tight')
