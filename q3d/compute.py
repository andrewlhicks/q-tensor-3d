from sympyplus import *

# functions that give the individual terms

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

def term_twist_var(q,p):
    """ Returns the variational derivative of the twist term """
    Q = QTensor(q)
    P = QTensor(p)

    return 2*c.q0*mixedp(Q,P) + 2*c.q0*mixedp(P,Q)

# compute function, which puts the PDE system together

def compute():
    from q3d.config import constants as c
    from q3d.config import settings

    # set up Qvector objects
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
    Drr = Param([Dr,r])
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

    Q0 = c.S0*(outerp(nu,nu) - (1.0/3.0)*eye(3))
    Pi = eye(3) - outerp(nu,nu)

    # Define 'energies' used to calculate the energy    
    energies = EnergyForm(Dqq,Dpp,Drr)
    energies.add_domain(GeneralForm(c.L1/2*termL1(Dq,Dq)+c.L2/2*termL2(Dq,Dq)+c.L3/2*termL3(Dq,Dq),Dqq,name='Elastic Energy'))
    if settings.pde.formulation == 'default':
        energies.add_domain(GeneralForm(2*c.q0*c.L1*mixedp(Q,Q),Dqq,name='Twist Energy'))
    elif settings.pde.formulation == 'lavrentovich':
        energies.add_domain(GeneralForm(2*c.q0*c.L1*mixedp(Q,Q) + 2*c.q0**2*c.L1*innerp(Q,Q),Dqq,name='Twist Energy'))
    energies.add_domain(GeneralForm((1/c.ep**2)*(1 - (c.A/2)*innerp(Q,Q) - (c.B/3)*trace(Q**3) + (c.C/4)*trace(Q**2)**2),Dqq,name='Bulk Energy'))
    energies.add_domain(GeneralForm(-f.dot(q),Dqq,name='Domain Forcing Energy'))
    energies.add_boundary(GeneralForm(c.W0/2*innerp(Q-Q0,Q-Q0),Dqq,name='Homeotropic Anchoring'))
    energies.add_boundary(GeneralForm(c.W1/2*innerp(Q-Pi*Q*Pi+c.S0/3*outerp(nu,nu),Q-Pi*Q*Pi+c.S0/3*outerp(nu,nu))+c.W2/4*(innerp(Q,Q)-2*c.S0**2/3)**2,Dqq,name='Planar-degenerate Anchoring'))
    energies.add_boundary(GeneralForm(-g.dot(q),Dqq,name='Boundary Forcing Energy'))

    energies_minmom = energies.copy()
    energies_minmom.add_domain(GeneralForm(1/(2*c.dt) * innerp(Q-QP,Q-QP),Dqq,name='Min Moments Energy'))

    # Construct PDE

    lhs_d = [ # lhs domain
        GeneralForm(c.L1*termL1(Dq,Dp)+c.L2*termL2(Dq,Dp)+c.L3*termL3(Dq,Dp),Dqq,Dpp,name='a_E(Q,P)'),
        GeneralForm(2*c.q0*c.L1*(mixedp(Q,P) + mixedp(P,Q)),Dqq,Dpp,name='a_T(Q,P)'),
        GeneralForm((1/c.ep**2)*(-c.A*innerp(Q,P) - c.B*innerp(Q**2,P) + c.C*innerp(Q,Q)*innerp(Q,P)),Dqq,Dpp,name='Dψ(Q):P')
        ]
    if settings.pde.formulation == 'lavrentovich':
        del lhs_d[1]
        lhs_d.append(GeneralForm(2*c.q0*c.L1*(mixedp(Q,P) + mixedp(P,Q)) + 4*c.q0**2*c.L1*innerp(Q,P),Dqq,Dpp,name='a_T(Q,P)'))
    rhs_d = [ # rhs domain
        GeneralForm(f.dot(p),Dpp,name='f(P)')
    ]
    lhs_b = [ # lhs boundary
        GeneralForm(c.W0*q.dot(p),Dqq,Dpp,name='W0(Q:P)'),
        GeneralForm(c.W1*innerp(Q-Pi*Q*Pi,P),Dqq,Dpp,name='W1(Q-ΠQΠ):P'),
        GeneralForm(c.W2*((innerp(Q,Q) - 2*c.S0**2/3)*innerp(Q,P)),Dqq,Dpp,name='(|Q|^2-2S0^2/3)(Q:P)')
    ]
    rhs_b = [ # rhs boundary
        GeneralForm(c.W0*innerp(Q0,P),Dpp,name='W0(Q0:P)'),
        GeneralForm(c.W1*innerp(-c.S0/3*outerp(nu,nu),P),Dpp,name='(-S0/3(nu⊗nu):P)'),
        GeneralForm(g.dot(p),Dpp,name='g(P)')
    ]

    if settings.pde.grad_desc:
        lhs_d.append(GeneralForm((1/c.dt)*(1/c.ep**2)*q.dot(p),Dqq,Dpp,name='(Ω)Q:P/dt'))
        lhs_b.append(GeneralForm((1/c.dt)*(c.W2)*q.dot(p),Dqq,Dpp,name='(Γ)Q:P/dt'))
        rhs_d.append(GeneralForm((1/c.dt)*(1/c.ep**2)*qp.dot(p),Dpp,name='(Ω)Q_p:P/dt'))
        rhs_b.append(GeneralForm((1/c.dt)*(c.W2)*qp.dot(p),Dpp,name='(Γ)Q_p:P/dt'))

    # assemble two PDEs
    pde_d = PDE(lhs_d,rhs_d,Dqq,Dpp,over='domain')
    pde_b = PDE(lhs_b,rhs_b,Dqq,Dpp,over='boundary')

    # apply newton's method
    pde_nm_d = pde_d.newtons_method(Dqnpqnp,Dpp,Dqq)
    pde_nm_b = pde_b.newtons_method(Dqnpqnp,Dpp,Dqq)

    # to make PDE system with positive definite lhs, remove non positive definite parts
    pde_pd_d = pde_nm_d.copy()
    pde_pd_d.rmv_lhs_form('∂(Dψ(Q):P)')
    pde_pd_d.rmv_lhs_form('∂(a_T(Q,P))')

    # Create relevant UFL strings

    out = {
        'pde_nm_d' : pde_nm_d.ufl,
        'pde_nm_b' : pde_nm_b.ufl,
        'energies' : energies.ufl_dict,
        'energies_minmom' : energies_minmom.ufl_dict,
        'pde' : [pde_d.ufl, pde_b.ufl],
        'pde_nm' : [pde_nm_d.ufl, pde_nm_b.ufl],
        'pde_pd' : [pde_pd_d.ufl, pde_nm_b.ufl],
    }

    return out

# END OF CODE
