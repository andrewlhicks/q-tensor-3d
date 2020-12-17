""" Here we write the code for computing the energy completely in UFL. """

from firedrake import *
from firedrakeplus import tensorfy
from settings import const
from sympyplus import levi_civita

def elastic(q):
	Q = tensorfy(q)
	summ = const.L1/2*Q[i,j].dx(k)*Q[i,j].dx(k) + const.L2/2*Q[i,j].dx(j)*Q[i,k].dx(k) + const.L3/2*Q[i,j].dx(k)*Q[i,k].dx(j)

	for ii in range(3):
		for jj in range(3):
			for kk in range(3):
				for ll in range(3):
					summ += const.L1/2*const.q0*levi_civita(ii,kk,ll)*Q[ll,jj].dx(kk)*Q[ii,jj]

	summ += const.L1/2*const.q0**2*inner(Q,Q)

	return summ

def bulk(q):
	Q = tensorfy(q)
	return 1/const.ep**2*(-const.A/2*inner(Q,Q) - const.B/3*inner(Q*Q,Q) + const.C/4*inner(Q,Q)**2)

def anchor_n(q,nu):
	Q = tensorfy(q)
	Q0 = const.S0*(outer(nu,nu)-1/3*Identity(3))
	return const.W0/2*inner(Q-Q0,Q-Q0)

def anchor_pd(q,nu):
	Q = tensorfy(q) + const.S0/3*Identity(3)
	Pi = Identity(3) - outer(nu,nu)
	return const.W1/2*inner(Q-Pi*Q*Pi,Q-Pi*Q*Pi) - const.W1/2*inner(Q,outer(nu,nu)) + const.W2/4*(inner(Q,Q)-const.S0**2)**2 # energy is "trace-free"
