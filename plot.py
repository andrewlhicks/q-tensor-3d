import matplotlib.pyplot as plt
import settings
import saves
import numpy as np

def time_vs_energy(times,energies,refinement_level='Not specified'):
	fig, ax = plt.subplots(figsize=(10,10))
	fig.suptitle(f'Energy over {settings.mesh.name} Mesh',fontsize=16)

	line_style = '--' if len(energies) > 50 else 'o--'

	plot_energies = saves.EnergyList(energies)

	plot_energies[0] = None
	for ii in range(len(plot_energies)):
		if ii != 0:
			plot_energies[ii] -= energies[ii-1]

	ax.plot(times,plot_energies,line_style)

	ax.set_xlabel('Time')
	ax.set_ylabel('Energy')
	# ax.set_yscale('log')
	ax.set_title(f'Refinement level: {refinement_level}')

	plt.savefig(f'{saves.current_directory}/energy/ref_{refinement_level}.png')
	plt.close()

def scatter_vs_poly(scatter,poly):
	fig, ax = plt.subplots(figsize=(10,10))
	fig.suptitle(f'Scatter vs. curve',fontsize=16)

	x = np.linspace(0,10,5)
	ax.plot(x,scatter,'ro')

	x = np.linspace(0,10,30)
	ax.plot(poly(x),'b')

	ax.set_xlabel('xi')
	ax.set_ylabel('E')

	plt.savefig(f'{saves.current_directory}/energy/poly.png')
	plt.close()
