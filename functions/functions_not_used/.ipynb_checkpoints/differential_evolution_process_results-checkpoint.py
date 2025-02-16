import numpy as np
import pandas as pd

def process_results(pop,energies):
    energies = np.array(energies)

    pop_up = [pop[0]]
    energies_up = [energies[0]]
    
    for i in range(1, len(pop)):
        new_energies = np.copy(energies_up[-1])
        new_pop = np.copy(pop_up[-1])
        pos = energies[i] < new_energies
        new_energies[pos] = energies[i][pos]
        new_pop[pos] = pop[i][pos]
        pop_up.append(new_pop)
        energies_up.append(new_energies)

    best_evo = pd.DataFrame(pop_up[-1])

    best_evo["penalty"]=energies_up[-1]

    best_evo_sorted = best_evo.sort_values("penalty")

    return best_evo_sorted