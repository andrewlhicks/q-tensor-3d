""" Here we write the code for computing the energy completely in UFL. """

from firedrake import *
from firedrakeplus import tensorfy
from settings import const

def elastic(q):
	Q = tensorfy(q)
	return const.L1/2*Q[i,j].dx(k)*Q[i,j].dx(k) + const.L2/2*Q[i,j].dx(j)*Q[i,k].dx(k) + const.L3/2*Q[i,j].dx(k)*Q[i,k].dx(j)

def bulk(q):
	Q = tensorfy(q)
	return 1/const.ep**2*(-const.A/2*inner(Q,Q) - const.B/3*inner(Q*Q,Q) + const.C/4*inner(Q,Q)**2)

def forcing_f(q,f):
	return dot(f,q)

def forcing_g(q,g):
	return dot(g,q)

def anchor_n(q,nu):
	Q = tensorfy(q)
	Q0 = const.S0*(outer(nu,nu)-1/3*Identity(3))
	return const.W0/2*inner(Q-Q0,Q-Q0)

def anchor_pd(q,nu):
	Q = tensorfy(q) + const.S0/3*Identity(3)
	Pi = Identity(3) - outer(nu,nu)
	return const.W1/2*inner(Q-Pi*Q*Pi,Q-Pi*Q*Pi) + const.W2/4*(inner(Q,Q)-const.S0**2)**2
