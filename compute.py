from sympyplus import *

import user_expressions.initial_q as initial_q
import user_expressions.manufac_q as manufac_q
import user_expressions.forcing_f as forcing_f
import user_expressions.forcing_g as forcing_g
import user_expressions.bdycond_s as bdycond_s

from config import settings
from config import constants as c

# Set up Qvector objects

nu = AbstractVector('nu')

q = QVector('q')
Dq = q.grad
Dqq = Param([Dq,q])
Q = q.tens

p = QVector('p')
Dp = p.grad
Dpp = Param([Dp,p])
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
Dqnpqnp = Param([Dqnp,qnp])
QNP = qnp.tens

f = QVector('f')
g = QVector('g')

c.S0 = 0.700005530940965

Q0 = c.S0*(outerp(nu,nu) - (1.0/3.0)*eye(3))
Pi = eye(3) - outerp(nu,nu)

# theta = atan2(x[1]-0.5,x[0]-0.5)
# N = Matrix([cos(theta),sin(theta),0])
# Q0 = c.S0*(outerp(N,N) - (1/3)*eye(3))

#################################################################

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

# def term_twist(q):
#     """ Returns the twist term from the energy """
#     Q = QTensor(q)

#     return c.q0*mixedp(Q,Q) + c.q0**2*innerp(Q,Q)

def term_twist_var(q,p):
    """ Returns the variational derivative of the twist term """
    Q = QTensor(q)
    P = QTensor(p)

    # return 2*c.q0*mixedp(Q,P) + 2*c.q0*mixedp(P,Q) + 4*c.q0**2*innerp(Q,P)
    return 2*c.q0*mixedp(Q,P) + 2*c.q0*mixedp(P,Q) # got rid of the linear 0-th derivative term

#########################################################

def compute():
    # Define 'energies' used to calculate the energy    
    energies = EnergyForm([Dq,q],[Dp,p],[Dr,r])
    energies.add_domain(GeneralForm(c.L1/2*termL1(Dq,Dq)+c.L2/2*termL2(Dq,Dq)+c.L3/2*termL3(Dq,Dq),[Dq,q],name='Elastic Energy'))
    energies.add_domain(GeneralForm(2*c.L1*c.q0*mixedp(Q,Q),[Dq,q],name='Twist Energy'))
    energies.add_domain(GeneralForm((1/c.ep**2)*(1 - (c.A/2)*innerp(Q,Q) - (c.B/3)*trace(Q**3) + (c.C/4)*trace(Q**2)**2),[Dq,q],name='Bulk Energy'))
    energies.add_domain(GeneralForm(-f.dot(q),[Dq,q],name='Domain Forcing Energy'))
    energies.add_boundary(GeneralForm(c.W0/2*innerp(Q-Q0,Q-Q0),[Dq,q],name='Homeotropic Anchoring'))
    energies.add_boundary(GeneralForm(c.W1/2*innerp(Q-Pi*Q*Pi+c.S0/3*outerp(nu,nu),Q-Pi*Q*Pi+c.S0/3*outerp(nu,nu))+c.W2/4*(innerp(Q,Q)-2*c.S0**2/3)**2,[Dq,q],name='Planar-degenerate Anchoring'))
    energies.add_boundary(GeneralForm(-g.dot(q),[Dq,q],name='Boundary Forcing Energy'))

    # Construct PDE

    lhs_d = [ # lhs domain
        GeneralForm((1/c.dt)*q.dot(p),[Dq,q],[Dp,p],name='Q:P/dt'),
        GeneralForm(c.L1*termL1(Dq,Dp)+c.L2*termL2(Dq,Dp)+c.L3*termL3(Dq,Dp),[Dq,q],[Dp,p],name='a_E(Q,P)'),
        GeneralForm(c.L1*term_twist_var(q,p),[Dq,q],[Dp,p],name='a_T(Q,P)'),
        GeneralForm((1/c.ep**2)*(-c.A*innerp(Q,P) - c.B*innerp(Q**2,P) + c.C*innerp(Q,Q)*innerp(Q,P)),[Dq,q],[Dp,p],name='Dψ(Q):P')
        ]
    rhs_d = [ # rhs domain
        GeneralForm((1/c.dt)*qp.dot(p),[Dp,p],name='Q_p:P/dt'),
        GeneralForm(f.dot(p),[Dp,p],name='f(P)')
    ]
    lhs_b = [ # lhs boundary
        GeneralForm(c.W0*q.dot(p),[Dq,q],[Dp,p],name='W0(Q:P)'),
        GeneralForm(c.W1*innerp(Q-Pi*Q*Pi,P),[Dq,q],[Dp,p],name='W1(Q-ΠQΠ):P'),
        GeneralForm(c.W2*((innerp(Q,Q) - 2*c.S0**2/3)*innerp(Q,P)),[Dq,q],[Dp,p],name='(|Q|^2-2S0^2/3)(Q:P)')
    ]
    rhs_b = [ # rhs boundary
        GeneralForm(c.W0*innerp(Q0,P),[Dp,p],name='W0(Q0:P)'),
        GeneralForm(c.W1*innerp(-c.S0/3*outerp(nu,nu),P),[Dp,p],name='(-S0/3(nu⊗nu):P)'),
        GeneralForm(g.dot(p),[Dp,p],name='g(P)')
    ]

    # assemble two PDEs
    pde_d = PDE(lhs_d,rhs_d,Dqq,Dpp,over='domain')
    pde_b = PDE(lhs_b,rhs_b,Dqq,Dpp,over='boundary')

    # apply newton's method
    pde_d = pde_d.newtons_method(Dqnpqnp,Dpp,Dqq)
    pde_b = pde_b.newtons_method(Dqnpqnp,Dpp,Dqq)

    # Create relevant UFL strings

    out = {
        'initial_q' : initial_q.out,
        'manufac_q' : manufac_q.out,
        'forcing_f' : forcing_f.out,
        'forcing_g' : forcing_g.out,
        'bdycond_s' : bdycond_s.out,
        'n_bf_O' : pde_d.uflfy()['lhs'],
        'n_bf_G' : pde_b.uflfy()['lhs'],
        'n_lf_O' : pde_d.uflfy()['rhs'],
        'n_lf_G' : pde_b.uflfy()['rhs'],
        'energies' : energies.ufl_dict
    }

    return out

# END OF CODE
