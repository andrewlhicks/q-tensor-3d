""" Treatment of the L1, L2, and L3 terms, since they are too complicated to write in-line. """

from sympyplus import *
from settings import const

nu = Matrix([x[0],x[1],x[2]])

def termL1(grad1,grad2):
    return innerp(grad1,grad2)

def termL2(grad1,grad2):
    term = 0

    for ii in range(5):
        for jj in range(5):
            term += (grad1.row(ii) * E[ii]).dot(E[jj] * grad2.row(jj).T)
    
    return term

def termL3(grad1,grad2):
    term = 0

    for ii in range(5):
        for jj in range(5):
            term += (grad1.row(ii) * E[jj]).dot(E[ii] * grad2.row(jj).T)
    
    return term

def term_twist(q):
    """ Returns the twist term from the energy """
    Q = QTensor(q)

    return const.q0*mixedp(Q,Q) + const.q0**2*innerp(Q,Q)

def term_twist_var(q,p):
    """ Returns the variational derivative of the twist term """
    Q = QTensor(q)
    P = QTensor(p)

    return const.q0*mixedp(Q,P) + const.q0*mixedp(P,Q) + 2*const.q0**2*innerp(Q,P)

####

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
                    term[ii,jj] -= const.q0*levi_civita(ii,ll,kk)*diff(Q[ll,jj],x[kk])

    term += const.q0**2*Q

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
                    term[ii,jj] += const.q0/2*nu[kk]*levi_civita(ii,ll,kk)*Q[ll,jj]

    return term

S0 = (const.B + sqrt(const.B**2 + 24.0*const.A*const.C))/(4.0*const.C)
Q0 = S0*(outerp(nu,nu) - (1.0/3.0)*eye(3))
Pi = eye(3) - outerp(nu,nu)

def tilde(Q):
    return Q + S0/3*eye(3)

###

def strongForm(G): # plugs G into the strong form PDE
    return const.L1*strongL1(G) + const.L2*strongL2(G) + const.L3*strongL3(G) + const.L1*strong_twist(G) + (1/const.ep**2)*(-const.A*G - const.B*G*G + const.C*innerp(G,G)*G)

def strongFormGamma(G):
    return const.L1*strongGammaL1(G) + const.L2*strongGammaL2(G) + const.L3*strongGammaL3(G) + const.L1*strong_twist_gamma(G) + const.W0*(G-Q0) + const.W1*(tilde(G)-Pi*tilde(G)*Pi) + const.W2*(innerp(tilde(G),tilde(G))-S0**2)*G