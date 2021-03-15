from sympyplus import *

import user_expressions.initial_q as initial_q
import user_expressions.manufac_q as manufac_q
import user_expressions.forcing_f as forcing_f
import user_expressions.forcing_g as forcing_g
import user_expressions.bdycond_s as bdycond_s
import user_expressions.bdycond_w as bdycond_w

def _set_settings_file(settings_file_path):
    import settings as st
    import constants as ct

    st._load_file(settings_file_path)
    global settings
    settings = st

    ct._load_file(settings.constants.file_path)
    global const
    const = ct.const
    const.dt = settings.timedata.time_step

    global S0, Q0, Pi
    S0 = (const.B + sqrt(const.B**2 + 24.0*const.A*const.C))/(4.0*const.C)
    Q0 = S0*(outerp(nu,nu) - (1.0/3.0)*eye(3))
    Pi = eye(3) - outerp(nu,nu)

# theta = atan2(x[1]-0.5,x[0]-0.5)
# N = Matrix([cos(theta),sin(theta),0])
# Q0 = const.S0*(outerp(N,N) - (1/3)*eye(3))

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

#     return const.q0*mixedp(Q,Q) + const.q0**2*innerp(Q,Q)

def term_twist_var(q,p):
    """ Returns the variational derivative of the twist term """
    Q = QTensor(q)
    P = QTensor(p)

    # return 2*const.q0*mixedp(Q,P) + 2*const.q0*mixedp(P,Q) + 4*const.q0**2*innerp(Q,P)
    return 2*const.q0*mixedp(Q,P) + 2*const.q0*mixedp(P,Q) # got rid of the linear 0-th derivative term

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

    # Got rid of this linear 0-th derivative term below:
    # term += const.q0**2*Q

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
                    term[ii,jj] += const.q0*nu[kk]*levi_civita(ii,ll,kk)*Q[ll,jj]

    return term

def tilde(Q):
    return Q + S0/3*eye(3)

###

def strong_F(Q_manu): # plugs Q_manu into the strong form PDE
    return const.L1*strongL1(Q_manu) + const.L2*strongL2(Q_manu) + const.L3*strongL3(Q_manu) + 4*const.L1*strong_twist(Q_manu) + (1/const.ep**2)*(-const.A*Q_manu - const.B*Q_manu*Q_manu + const.C*innerp(Q_manu,Q_manu)*Q_manu)

def strong_G(Q_manu):
    return const.L1*strongGammaL1(Q_manu) + const.L2*strongGammaL2(Q_manu) + const.L3*strongGammaL3(Q_manu) + 2*const.L1*strong_twist_gamma(Q_manu) + const.W0*(Q_manu-Q0) + const.W1*(tilde(Q_manu)-Pi*tilde(Q_manu)*Pi-trace(tilde(Q_manu)-Pi*tilde(Q_manu)*Pi)/3*eye(3)) + const.W2*(innerp(tilde(Q_manu),tilde(Q_manu))-S0**2)*Q_manu

#########################################################

def compute():
    # Energies

    energyElastic = GeneralForm(const.L1/2*termL1(Dq,Dq)+const.L2/2*termL2(Dq,Dq)+const.L3/2*termL3(Dq,Dq),[Dq,q],name='energyElastic')
    energyBulkC = GeneralForm((1/const.ep**2)*(((const.L0 - const.A)/2)*innerp(Q,Q) - (const.B/3)*trace(Q**3) + (const.C/4)*trace(Q**2)**2),[Dq,q])
    energyBulkE = GeneralForm((1/const.ep**2)*(const.L0/2)*trace(Q**2),[Dq,q])

    # Bilinear forms

    bfTimeStep = GeneralForm((1/const.dt)*q.dot(p),[Dq,q],[Dp,p],name='bfTimeStep')
    bfElastic = variationalDerivative(energyElastic,[Dq,q],[Dp,p],name='bfElastic')
    bfBulkC = variationalDerivative(energyBulkC,[Dq,q],[Dp,p],name='bfBulkC')
    bfBulkE = variationalDerivative(energyBulkE,[Dq,q],[Dp,p],name='bfBulkE')

    bfTwist = GeneralForm(const.L1*term_twist_var(q,p),[Dq,q],[Dp,p],name='bfTwist')

    bfNAnchor = GeneralForm(const.W0*q.dot(p),[Dq,q],[Dp,p],name='bfNAnchor')
    bfPDAnchor1 = GeneralForm(const.W1*innerp(Q-Pi*Q*Pi,P),[Dq,q],[Dp,p],name='bfPDAnchor1') # This was WRONG and nonlinear. Why didn't the code catch this?
    bfPDAnchor2 = GeneralForm(const.W2*((innerp(Q,Q) - 2*S0**2/3)*innerp(Q,P)),[Dq,q],[Dp,p],name='bfPDAnchor2')

    # Linear forms

    lfTimeStep = GeneralForm(bfTimeStep([Dqp,qp],[Dp,p]),[Dp,p],name='lfTimeStep')
    lfBulkE = GeneralForm(bfBulkE([Dqp,qp],[Dp,p]),[Dp,p],name='lfBulkE')

    lfNAnchor = GeneralForm(const.W0*innerp(Q0,P),[Dp,p],name='lfNAnchor')
    lfPDAnchor1 = GeneralForm(const.W1*innerp(-S0/3*outerp(nu,nu),P),[Dp,p],name='lfPDAnchor1')

    lf_ForcingF = GeneralForm(f.dot(p),[Dp,p],name='lf_ForcingF')
    lf_ForcingG = GeneralForm(g.dot(p),[Dp,p],name='lf_ForcingG')

    # Assemble LHS, RHS

    bilinearDomain = lhsForm([Dq,q],[Dp,p],forms=[bfTimeStep,bfElastic,bfBulkC,bfTwist])
    bilinearBoundary = lhsForm([Dq,q],[Dp,p],forms=[bfNAnchor,bfPDAnchor1,bfPDAnchor2])

    linearDomain = rhsForm([Dp,p],forms=[lfTimeStep,lfBulkE,lf_ForcingF])
    linearBoundary = rhsForm([Dp,p],forms=[lfNAnchor,lfPDAnchor1,lf_ForcingG])

    # Newton's method

    newt_bilinearDomain, newt_linearDomain = newtonsMethod(bilinearDomain,linearDomain,[Dqnp,qnp],[Dq,q],[Dp,p])
    newt_bilinearBoundary, newt_linearBoundary = newtonsMethod(bilinearBoundary,linearBoundary,[Dqnp,qnp],[Dq,q],[Dp,p])

    # Create relevant UFL strings

    class out:
        initial_q = initial_q.out
        manufac_q = manufac_q.out
        forcing_f = forcing_f.out
        forcing_g = forcing_g.out
        bdycond_s = bdycond_s.out
        bdycond_w = bdycond_w.out

        n_bf_O = newt_bilinearDomain()
        n_bf_G = newt_bilinearBoundary()
        n_lf_O = newt_linearDomain()
        n_lf_G = newt_linearBoundary()

    return out

# END OF CODE