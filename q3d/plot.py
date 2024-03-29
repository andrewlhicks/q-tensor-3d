import matplotlib.pyplot as plt
import numpy as np
import os
import q3d.saves as saves

def usage():
	usage = """usage: python plot.py [-o] (-l | -r) <save-name>
  l: file in './saves'
  r: file in './saves-remote'
  o: open file"""
	print(usage)

def time_vs_energy(times,energies,refinement_level='Not specified',open_file=False):
	from q3d.config import settings
	import matplotlib

	if not open_file:
		matplotlib.use('Agg')

	try:
		fig, (ax1, ax2) = plt.subplots(2,1,figsize=(10,10))
		fig.suptitle(f'{settings.mesh.name} Mesh, Ref. {refinement_level}',fontsize=16)

		line_style = '--' if len(energies) > 50 else 'o--'

		plot_energies = saves.EnergyList(energies)

		plot_energies[0] = None
		for ii in range(len(plot_energies)):
			if ii != 0:
				plot_energies[ii] -= energies[ii-1]

		# Truncate to last 1000

		end_value = times[-1]
		N = 100000
		start_value = np.maximum(end_value-N*settings.time.step,0)

		ax1.plot(times[-N:],energies[-N:],line_style)
		ax2.plot(times[-N:],plot_energies[-N:],line_style,color='red')

		ax1.set_xlabel('Time')
		ax1.set_xlim(start_value,end_value)
		ax1.set_ylabel('Energy')
		ax1.set_title('Total energy')
		ax1.grid()

		ax2.set_xlabel('Time')
		ax2.set_xlim(start_value,end_value)
		ax2.set_ylabel('Energy')
		ax2.set_title('Change in energy')
		ax2.set_yscale('symlog')
		ax2.grid()

		file_path = f'{saves.SavePath}/energy/ref_{refinement_level}.png'
		
		plt.savefig(file_path)

		if open_file:
			plt.show()
		
		plt.close()

		return file_path
	except Exception as err: # temporary replacement for TclError until I figure out why that doesn't work anymore
		print(err)

def scatter_vs_poly(scatter,poly):
	fig, ax = plt.subplots(figsize=(10,10))
	fig.suptitle(f'Scatter vs. curve',fontsize=16)

	x = np.linspace(0,10,5)
	ax.plot(x,scatter,'ro')

	x = np.linspace(0,10,30)
	ax.plot(poly(x),'b')

	ax.set_xlabel('xi')
	ax.set_ylabel('E')

	plt.savefig(f'{saves.SavePath}/energy/poly.png')
	plt.close()

def main():
	import getopt
	import sys
	import q3d.config as config

	try:
		opts, args = getopt.getopt(sys.argv[1:], 'l:r:o', ['help'])
	except getopt.GetoptError as err:
		print(err)  # will print something like "option -a not recognized"
		sys.exit()

	open_file = False

	for o, a in opts:
		if o in ('--help'):
			usage()
			sys.exit()
		elif o in ('-l'):
			remote = ''
			save_name = a
			saves.initialize('',save_name)
			config.initialize(f'saves/{save_name}/settings.yml') # So plot.py has mesh name
		elif o in ('-r'):
			remote = 'remote '
			save_name = a
			saves.initialize('',save_name,remote=True)
			config.initialize(f'saves-remote/{save_name}/settings.yml') # So plot.py has mesh name
		elif o in ('-o'):
			open_file = True
		else:
			assert False, "unhandled option"

	try:
		save_name
	except NameError:
		print("Must choose local (-l) or remote (-r).")
		sys.exit()

	times, energies = saves.load_energies()
	file_path = time_vs_energy(times,energies,open_file=open_file)
	print(f"Done plotting to {remote}save {save_name}.")

if __name__ == '__main__':
	main()
