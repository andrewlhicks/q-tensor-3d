# We must re-import "settings" here, because Python is weird like that and it won't work otherwise

import sympy as sp
import sympyplus as spp
from settings import *

def boundary():
    n = bound_cond
    G = spp.outerp(n,n) - (1.0/3.0) * sp.eye(3)
    return spp.uflfy(spp.vectorfy(G))

def bilinear():
    a = (sp.sqrt(3.0)-3.0)/6.0
    b = (sp.sqrt(3.0)+3.0)/6.0
    c = -sp.sqrt(3.0)/3.0
    d = sp.sqrt(2.0)/2.0

    E = [sp.diag(a,b,c), sp.diag(b,a,c), sp.Matrix([[0,d,0],[d,0,0],[0,0,0]]), sp.Matrix([[0,0,d],[0,0,0],[d,0,0]]), sp.Matrix([[0,0,0],[0,0,d],[0,d,0]])]
    
    vrs = [x,y,z]
    
    q0, q1, q2, q3, q4 = sp.symbols('q[0:5]')

    d0q0, d1q0, d2q0, \
    d0q1, d1q1, d2q1, \
    d0q2, d1q2, d2q2, \
    d0q3, d1q3, d2q3, \
    d0q4, d1q4, d2q4 = sp.symbols('q[0:5].dx((0:3))')

    p0, p1, p2, p3, p4 = sp.symbols('p[0:5]')

    d0p0, d1p0, d2p0, \
    d0p1, d1p1, d2p1, \
    d0p2, d1p2, d2p2, \
    d0p3, d1p3, d2p3, \
    d0p4, d1p4, d2p4 = sp.symbols('p[0:5].dx((0:3))')
    
    qv = sp.Matrix([q0,q1,q2,q3,q4])
    pv = sp.Matrix([p0,p1,p2,p3,p4])
    
    grad_q = sp.Matrix([[d0q0,d1q0,d2q0],
                         [d0q1,d1q1,d2q1],
                         [d0q2,d1q2,d2q2],
                         [d0q3,d1q3,d2q3],
                         [d0q4,d1q4,d2q4]])
                         
    grad_p = sp.Matrix([[d0p0,d1p0,d2p0],
                         [d0p1,d1p1,d2p1],
                         [d0p2,d1p2,d2p2],
                         [d0p3,d1p3,d2p3],
                         [d0p4,d1p4,d2p4]])

    term_L1 = spp.innerp(grad_q,grad_p)

    term_L2 = 0

    for ii in range(5):
        for jj in range(5):
            term_L2 += (grad_q.row(ii) * E[ii]).dot(E[jj] * grad_p.row(jj).T)

    term_L3 = 0

    for ii in range(5):
        for jj in range(5):
            term_L3 += (grad_q.row(ii) * E[jj]).dot(E[ii] * grad_p.row(jj).T)
    
    return spp.uflfy((1/dt) * qv.dot(pv) + L1*term_L1 + L2*term_L2 + L3*term_L3 + (L0/ep)*qv.dot(pv))

def linear():
    a = (sp.sqrt(3.0)-3.0)/6.0
    b = (sp.sqrt(3.0)+3.0)/6.0
    c = -sp.sqrt(3.0)/3.0
    d = sp.sqrt(2.0)/2.0

    E = [sp.diag(a,b,c), sp.diag(b,a,c), sp.Matrix([[0,d,0],[d,0,0],[0,0,0]]), sp.Matrix([[0,0,d],[0,0,0],[d,0,0]]), sp.Matrix([[0,0,0],[0,0,d],[0,d,0]])]
    
    vrs = [x,y,z]
    
    qp0, qp1, qp2, qp3, qp4 = sp.symbols('q_prev[0:5]')
    qpv = sp.Matrix([qp0,qp1,qp2,qp3,qp4])
    
    p0, p1, p2, p3, p4 = sp.symbols('p[0:5]')
    pv = sp.Matrix([p0,p1,p2,p3,p4])
    
    if manufactured == 0:
        fv = sp.Matrix([0,0,0,0,0])
    elif manufactured == 1:
        n = bound_cond
        Gt = spp.outerp(n,n) - (1.0/3.0) * sp.eye(3)
        gv = spp.vectorfy(Gt)
        
        strong_L1 = sp.zeros(3,3)
        
        for ii in range(3):
            for jj in range(3):
                for kk in range(3):
                    strong_L1[ii,jj] -= sp.diff(Gt[ii,jj],vrs[kk],vrs[kk])
        
        strong_L2 = sp.zeros(3,3)
        
        for ii in range(3):
            for jj in range(3):
                for kk in range(3):
                    strong_L2[ii,jj] -= sp.diff(Gt[ii,kk],vrs[jj],vrs[kk])
        
        strong_L3 = sp.zeros(3,3)
        
        for ii in range(3):
            for jj in range(3):
                for kk in range(3):
                    strong_L3[ii,jj] -= sp.diff(Gt[ii,kk],vrs[kk],vrs[jj])
        
        strong_A = Gt
        
        strong_B = Gt*Gt
        
        strong_C = spp.innerp(Gt,Gt)*Gt
        
        Ft = L1*strong_L1 + L2*strong_L2 + L3*strong_L3 + (1/ep)*(-A*strong_A - B*strong_B + C*strong_C)
        
        fv = spp.vectorfy(Ft)
    else:
        raise ValueError("Variable 'manufactured' must be 0 or 1.")
        
    term_A = qpv.dot(pv)

    Qpt = sp.zeros(3,3)
    Pt = sp.zeros(3,3)

    for ii in range(5):
        Qpt += qpv[ii]*E[ii]
        Pt += pv[ii]*E[ii]

    term_B = spp.innerp(Qpt*Qpt,Pt)

    term_C = qpv.dot(qpv) * qpv.dot(pv)
    
    return spp.uflfy((1/dt) * qpv.dot(pv) + (1/ep) * ( (A + L0) * term_A + B * term_B - C * term_C ) + fv.dot(pv))
    