import getopt
import matplotlib.pyplot as plt
import numpy as np
import saves
import settings
import sys
import yaml

def time_vs_energy(times,energies,refinement_level='Not specified'):
	fig, (ax1, ax2) = plt.subplots(2,1,figsize=(10,10))
	fig.suptitle(f'{settings.mesh.name} Mesh, Ref. {refinement_level}',fontsize=16)

	line_style = '--' if len(energies) > 50 else 'o--'

	plot_energies = saves.EnergyList(energies)

	plot_energies[0] = None
	for ii in range(len(plot_energies)):
		if ii != 0:
			plot_energies[ii] -= energies[ii-1]

	ax1.plot(times,energies,line_style)
	ax2.plot(times,plot_energies,line_style,color='red')

	ax1.set_xlabel('Time')
	ax1.set_xlim(0,len(times))
	ax1.set_ylabel('Energy')
	ax1.set_title('Total energy')

	ax2.set_xlabel('Time')
	ax2.set_xlim(0,len(times))
	ax2.set_ylabel('Energy')
	ax2.set_title('Change in energy')

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

def main():
	def usage():
		print("usage: python plot.py <save-name>")

	try:
		opts, args = getopt.getopt(sys.argv[1:], 'hl:r:', ['help,local=,remote='])
	except getopt.GetoptError as err:
		print(err)  # will print something like "option -a not recognized"
		sys.exit()

	for o, a in opts:
		if o in ('-h','--help'):
			usage()
			sys.exit()
		elif o in ('-l','--local'):
			remote = ''
			save_name = a
			saves.initialize(None,save_name)
			settings._load_file(f'saves/{save_name}/settings.yml')
		elif o in ('-r','--remote'):
			remote = 'remote '
			save_name = a
			saves.initialize(None,save_name,remote=True)
			settings._load_file(f'saves-remote/{save_name}/settings.yml')
		else:
			assert False, "unhandled option"

	times, energies = saves.load_energies()
	time_vs_energy(times,energies)
	print(f"Done plotting to {remote}save {save_name}.")

if __name__ == '__main__':
	main()
