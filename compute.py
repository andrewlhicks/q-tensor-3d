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
    # Classes for direct computation of energy and its derivatives
    # class energies:
    #     domain = [
    #         GeneralForm(c.L1/2*termL1(Dq,Dq)+c.L2/2*termL2(Dq,Dq)+c.L3/2*termL3(Dq,Dq),[Dq,q],name='Elastic Energy'),
    #         GeneralForm(2*c.L1*c.q0*mixedp(Q,Q),[Dq,q],name='Twist Energy'),
    #         # Standard double well
    #         # GeneralForm((1/c.ep**2)*(- (c.A/2)*innerp(Q,Q) - (c.B/3)*trace(Q**3) + (c.C/4)*trace(Q**2)**2),[Dq,q],name='Bulk Energy'),
    #         # Alternative double well
    #         GeneralForm((1/c.ep**2)*(1 - (c.A/2)*innerp(Q,Q) - (c.B/3)*trace(Q**3) + (c.C/4)*trace(Q**2)**2),[Dq,q],name='Bulk Energy'),
    #         GeneralForm(-f.dot(q),[Dq,q],name='Domain Forcing Energy')
    #     ]
    #     boundary = [
    #         GeneralForm(c.W0/2*innerp(Q-Q0,Q-Q0),[Dq,q],name='Homeotropic Anchoring'),
    #         GeneralForm(c.W1/2*innerp(Q-Pi*Q*Pi+c.S0/3*outerp(nu,nu),Q-Pi*Q*Pi+c.S0/3*outerp(nu,nu))+c.W2/4*(innerp(Q,Q)-2*c.S0**2/3)**2,[Dq,q],name='Planar-degenerate Anchoring'),
    #         GeneralForm(-g.dot(q),[Dq,q],name='Boundary Forcing Energy')
    #     ]
    # class energies_der:
    #     domain = [variationalDerivative(general_form,[Dq,q],[Dp,p]) for general_form in energies.domain]
    #     boundary = [variationalDerivative(general_form,[Dq,q],[Dp,p]) for general_form in energies.boundary]
    # class energies_der_der:
    #     domain = [secondVariationalDerivative(general_form,[Dq,q],[Dr,r],[Dp,p]) for general_form in energies_der.domain]
    #     boundary = [secondVariationalDerivative(general_form,[Dq,q],[Dr,r],[Dp,p]) for general_form in energies_der.boundary]
    
    energies = EnergyForm([Dq,q],[Dp,p],[Dr,r])
    energies.add_domain(GeneralForm(c.L1/2*termL1(Dq,Dq)+c.L2/2*termL2(Dq,Dq)+c.L3/2*termL3(Dq,Dq),[Dq,q],name='Elastic Energy'))
    energies.add_domain(GeneralForm(2*c.L1*c.q0*mixedp(Q,Q),[Dq,q],name='Twist Energy'))
    energies.add_domain(GeneralForm((1/c.ep**2)*(1 - (c.A/2)*innerp(Q,Q) - (c.B/3)*trace(Q**3) + (c.C/4)*trace(Q**2)**2),[Dq,q],name='Bulk Energy'))
    energies.add_domain(GeneralForm(-f.dot(q),[Dq,q],name='Domain Forcing Energy'))
    energies.add_boundary(GeneralForm(c.W0/2*innerp(Q-Q0,Q-Q0),[Dq,q],name='Homeotropic Anchoring'))
    energies.add_boundary(GeneralForm(c.W1/2*innerp(Q-Pi*Q*Pi+c.S0/3*outerp(nu,nu),Q-Pi*Q*Pi+c.S0/3*outerp(nu,nu))+c.W2/4*(innerp(Q,Q)-2*c.S0**2/3)**2,[Dq,q],name='Planar-degenerate Anchoring'))
    energies.add_boundary(GeneralForm(-g.dot(q),[Dq,q],name='Boundary Forcing Energy'))

    # Energies

    energyElastic = GeneralForm(c.L1/2*termL1(Dq,Dq)+c.L2/2*termL2(Dq,Dq)+c.L3/2*termL3(Dq,Dq),[Dq,q],name='energyElastic')
    # Standard double well
    # energyBulkC = GeneralForm((1/c.ep**2)*(((c.L0 - c.A)/2)*innerp(Q,Q) - (c.B/3)*trace(Q**3) + (c.C/4)*trace(Q**2)**2),[Dq,q])
    # energyBulkE = GeneralForm((1/c.ep**2)*(c.L0/2)*trace(Q**2),[Dq,q])
    # Alternative double well
    # energyBulkC = GeneralForm((1/c.ep**2)*(1 + (c.L0 - c.A)/2*trace(Q**2)),[Dq,q])
    # energyBulkE = GeneralForm((1/c.ep**2)*((c.L0/2)*innerp(Q,Q) + (c.B/3)*trace(Q**3) - (c.C/4)*trace(Q**2)**2),[Dq,q])
    # Fully implicit double well
    energyBulk = GeneralForm((1/c.ep**2)*(-(c.A/2)*innerp(Q,Q) - (c.B/3)*trace(Q**3) + (c.C/4)*trace(Q**2)**2),[Dq,q])

    # Bilinear forms

    bfTimeStep = GeneralForm((1/c.dt)*q.dot(p),[Dq,q],[Dp,p],name='bfTimeStep')
    bfElastic = variationalDerivative(energyElastic,[Dq,q],[Dp,p],name='bfElastic')
    # bfBulkC = variationalDerivative(energyBulkC,[Dq,q],[Dp,p],name='bfBulkC')
    # bfBulkE = variationalDerivative(energyBulkE,[Dq,q],[Dp,p],name='bfBulkE')
    bfBulk = variationalDerivative(energyBulk,[Dq,q],[Dp,p],name='bfBulk')

    bfTwist = GeneralForm(c.L1*term_twist_var(q,p),[Dq,q],[Dp,p],name='bfTwist')

    bfNAnchor = GeneralForm(c.W0*q.dot(p),[Dq,q],[Dp,p],name='bfNAnchor')
    bfPDAnchor1 = GeneralForm(c.W1*innerp(Q-Pi*Q*Pi,P),[Dq,q],[Dp,p],name='bfPDAnchor1') # This was WRONG and nonlinear. Why didn't the code catch this?
    bfPDAnchor2 = GeneralForm(c.W2*((innerp(Q,Q) - 2*c.S0**2/3)*innerp(Q,P)),[Dq,q],[Dp,p],name='bfPDAnchor2')

    # Linear forms

    lfTimeStep = GeneralForm(bfTimeStep([Dqp,qp],[Dp,p]),[Dp,p],name='lfTimeStep')
    # lfBulkE = GeneralForm(bfBulkE([Dqp,qp],[Dp,p]),[Dp,p],name='lfBulkE')

    lfNAnchor = GeneralForm(c.W0*innerp(Q0,P),[Dp,p],name='lfNAnchor')
    lfPDAnchor1 = GeneralForm(c.W1*innerp(-c.S0/3*outerp(nu,nu),P),[Dp,p],name='lfPDAnchor1')

    lf_ForcingF = GeneralForm(f.dot(p),[Dp,p],name='lf_ForcingF')
    lf_ForcingG = GeneralForm(g.dot(p),[Dp,p],name='lf_ForcingG')

    lf_TimeStep2 = GeneralForm( (c.beta)*(1/c.dt) * (qp.dot(p) - qpp.dot(p)) , [Dp,p], name='lf_TimeStep2')

    # Assemble LHS, RHS

    bilinearDomain = lhsForm([Dq,q],[Dp,p],name='Bilinear Domain',forms=[bfTimeStep,bfElastic,bfBulk,bfTwist])
    bilinearBoundary = lhsForm([Dq,q],[Dp,p],name='Bilinear Boundary',forms=[bfNAnchor,bfPDAnchor1,bfPDAnchor2])

    linearDomain = rhsForm([Dp,p],name='Linear Domain',forms=[lfTimeStep,lf_ForcingF])
    linearBoundary = rhsForm([Dp,p],name='Linear Boundary',forms=[lfNAnchor,lfPDAnchor1,lf_ForcingG])

    # If gradient descent is modified to the Heavy Ball method

    if settings.solver.grad_desc == 'heavyball':
        linearDomain.add_form(lf_TimeStep2)

    # Newton's method

    newt_bilinearDomain, newt_linearDomain = newtonsMethod(bilinearDomain,linearDomain,[Dqnp,qnp],[Dq,q],[Dp,p])
    newt_bilinearBoundary, newt_linearBoundary = newtonsMethod(bilinearBoundary,linearBoundary,[Dqnp,qnp],[Dq,q],[Dp,p])

    # Create relevant UFL strings

    out = {
        'initial_q' : initial_q.out,
        'manufac_q' : manufac_q.out,
        'forcing_f' : forcing_f.out,
        'forcing_g' : forcing_g.out,
        'bdycond_s' : bdycond_s.out,
        'n_bf_O' : newt_bilinearDomain(),
        'n_bf_G' : newt_bilinearBoundary(),
        'n_lf_O' : newt_linearDomain(),
        'n_lf_G' : newt_linearBoundary(),
        'energies' : energies.ufl_dict
    }

    return out

# END OF CODE
