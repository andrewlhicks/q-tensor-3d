import matplotlib.pyplot as plt
import settings
import saves

def time_vs_energy(times,energies,refinement_level='Not specified'):
	fig, ax = plt.subplots(figsize=(10,10))
	fig.suptitle(f'Energy over {settings.mesh.name} Mesh',fontsize=16)

	line_style = '--' if len(energies) > 50 else 'o--'

	ax.plot(times,energies,line_style)

	ax.set_xlabel('Time')
	ax.set_ylabel('Energy')
	ax.set_title(f'Refinement level: {refinement_level}')

	plt.savefig(f'{saves.current_directory}/energy/ref_{refinement_level}.png')
