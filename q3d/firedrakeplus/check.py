import q3d.printoff as pr

def check_energy_decrease(energies,current_time):
    """ Checks for energy decreasee in latest energy iteration. """
    if len(energies) < 2: return # Temporary escape clause until I can fix the EnergyList class
    change_in_energy = energies[-1]-energies[-2]
    energy_index = len(energies) # Will need to be fixed once I fix the EnergyList class
    threshold = 1.0e-12
    if change_in_energy > threshold:
        pr.warning(f'ΔE=+{change_in_energy:.05e} @t={current_time:.2f} @k={energy_index}')