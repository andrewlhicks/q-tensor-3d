import matplotlib.pyplot as plt

def time_vs_energy(ax,time,energy,manu_energy=None,refinement_level='Not specified'):
	ax.plot(time,energy,'--o')
	
	if manu_energy is not None:
		ax.hlines(manu_energy,time[0],time[-1])
	for i in range(1,len(energy)):
		ax.text(time[i],energy[i],energy[i],rotation=45)

	ax.set_xlabel('Time')
	ax.set_ylabel('Energy')
	ax.set_title(f'Refinement level: {refinement_level}')