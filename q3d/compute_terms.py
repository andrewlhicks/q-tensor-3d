""" This file's only purpose is to compute the strong form on Omega and on
Gamma. But since it defines its constants and its S0, Q0, and Pi
separately from compute.py, this file must eventually be eliminated. It is
an appendage that is only causing harm to the code as a whole. """

from sympyplus import *

from q3d.config import constants as c

# Set up Qvector objects

nu = AbstractVector('nu')

q = QVector('q')
Dq = q.grad
Q = q.tens

p = QVector('p')
Dp = p.grad
P = p.tens

r = QVector('r')
Dr = r.grad
R = r.tens

qp = QVector('q_prev')
Dqp = qp.grad
QP = qp.tens

qpp = QVector('q_prev_prev')
Dqpp = qpp.grad
QPP = qpp.tens

qnp = QVector('q_newt_prev')
Dqnp = qnp.grad
QNP = qnp.tens

f = QVector('f')
g = QVector('g')

Q0 = c.S0*(outerp(nu,nu) - (1.0/3.0)*eye(3))
Pi = eye(3) - outerp(nu,nu)

# Strong terms

def strongL1(Q):
    term = zeros(3,3)
    
    for ii in range(3):
            for jj in range(3):
                for kk in range(3):
                    term[ii,jj] -= diff(Q[ii,jj],x[kk],x[kk])
    
    return term

def strongL2(Q):
    term = zeros(3,3)
        
    for ii in range(3):
        for jj in range(3):
            for kk in range(3):
                term[ii,jj] -= diff(Q[ii,kk],x[jj],x[kk])
    
    return term

def strongL3(Q):  
    term = zeros(3,3)
        
    for ii in range(3):
        for jj in range(3):
            for kk in range(3):
                term[ii,jj] -= diff(Q[ii,kk],x[kk],x[jj])
    
    return term

def strong_twist(Q):
    term = zeros(3,3)

    for ii in range(3):
        for jj in range(3):
            for kk in range(3):
                for ll in range(3):
                    term[ii,jj] -= c.q0*levi_civita(ii,ll,kk)*diff(Q[ll,jj],x[kk])

    # Got rid of this linear 0-th derivative term below:
    # term += c.q0**2*Q

    return term

def strongGammaL1(Q):
    term = zeros(3,3)
    
    for ii in range(3):
            for jj in range(3):
                for kk in range(3):
                    term[ii,jj] += nu[kk]*diff(Q[ii,jj],x[kk])
    
    return term

def strongGammaL2(Q):
    term = zeros(3,3)
        
    for ii in range(3):
        for jj in range(3):
            for kk in range(3):
                term[ii,jj] += nu[jj]*diff(Q[ii,kk],x[kk])
    
    return term

def strongGammaL3(Q):  
    term = zeros(3,3)
        
    for ii in range(3):
        for jj in range(3):
            for kk in range(3):
                term[ii,jj] += nu[kk]*diff(Q[ii,kk],x[jj])
    
    return term

def strong_twist_gamma(Q):
    term = zeros(3,3)

    for ii in range(3):
        for jj in range(3):
            for kk in range(3):
                for ll in range(3):
                    term[ii,jj] += c.q0*nu[kk]*levi_civita(ii,ll,kk)*Q[ll,jj]

    return term

def tilde(Q):
    return Q + c.S0/3*eye(3)

###

def strong_F(Q_manu): # plugs Q_manu into the strong form PDE
    return c.L1*strongL1(Q_manu) + c.L2*strongL2(Q_manu) + c.L3*strongL3(Q_manu) + 4*c.L1*strong_twist(Q_manu) + (1/c.ep**2)*(-c.A*Q_manu - c.B*Q_manu*Q_manu + c.C*innerp(Q_manu,Q_manu)*Q_manu)

def strong_G(Q_manu):
    return c.L1*strongGammaL1(Q_manu) + c.L2*strongGammaL2(Q_manu) + c.L3*strongGammaL3(Q_manu) + 2*c.L1*strong_twist_gamma(Q_manu) + c.W0*(Q_manu-Q0) + c.W1*(tilde(Q_manu)-Pi*tilde(Q_manu)*Pi-trace(tilde(Q_manu)-Pi*tilde(Q_manu)*Pi)/3*eye(3)) + c.W2*(innerp(tilde(Q_manu),tilde(Q_manu))-c.S0**2)*Q_manu