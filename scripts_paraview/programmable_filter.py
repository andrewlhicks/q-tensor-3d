# set OUTWARD to True in order to make the normal vector be the "outward" normal (x,y,z) (for the SHELL)
# set OUTWARD to False in order to make the normal vector be (0,0,1) (for the SLAB)
OUTWARD = True

import numpy as np

# Get the first input
input = inputs[0]

# form Q-tensor
N = input.PointData['Q_00'].shape[0]
Q_tensor = np.zeros((N,3,3))

Q_tensor[:,0,0] = input.PointData['Q_00']
Q_tensor[:,1,0] = input.PointData['Q_10']
Q_tensor[:,2,0] = input.PointData['Q_20']
Q_tensor[:,2,1] = input.PointData['Q_21']
Q_tensor[:,2,2] = input.PointData['Q_22']

Q_tensor[:,0,1] = Q_tensor[:,1,0]
Q_tensor[:,0,2] = Q_tensor[:,2,0]
Q_tensor[:,1,2] = Q_tensor[:,2,1]

Q_tensor[:,1,1] = -(Q_tensor[:,0,0] + Q_tensor[:,2,2])

# obtain eigenvalues, eigenvectors
evals = np.zeros((N,3))
evecs = np.zeros((N,3,3))

for n in range(N):
    eval_vec_n, evec_mat_n = np.linalg.eigh(Q_tensor[n,:,:]) # extract eigenvalues, eigenvectors
    evals[n,:] = eval_vec_n[::-1] # place the largest eigenvalue first
    evecs[n,:,:] = evec_mat_n[:,::-1] # place the largest eigenvector first

# [output.PointData.append(evals[:,j], f'eval_{j}') for j in range(3)]
# [output.PointData.append(evecs[:,:,j], f'evec_{j}') for j in range(3)]

output.PointData.append(evecs[:,0], 'n')

# outward
outward = input.PointData['outward']

if OUTWARD:
    normal = outward
else:
    normal = np.zeros((N,3))
    for n in range(N):
        normal[n,:] = np.array([0,0,1])

output.PointData.append(normal, 'ν')

# compute 'magnitude'
mag = np.zeros((N,))
for n in range(N):
    mag[n] = abs(np.dot(normal[n,:],evecs[n,:,0]))

output.PointData.append(mag, 'abs(n•ν|)')

# compute norm of Q
norm = np.sqrt(sum(Q_tensor[:,i,j]**2 for i, j in zip(range(3),range(3))))

output.PointData.append(norm, '|Q|')