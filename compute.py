from sympyplus import *
from compute_terms import *
from settings import const, options

# Initial guess and manufactured solution

# n = Matrix([cos(x[0]+x[1]+x[2]),sin(x[0]+x[1]+x[2]),0])
# G = outerp(n,n) - (1.0/3.0) * eye(3)
X = Matrix([[cos(x[0]),sin(x[1]),cos(x[2])],
			[sin(x[1]),cos(x[1]),sin(x[2])],
			[cos(x[2]),sin(x[2]),sin(x[0])]])
M = X - trace(X)/3*eye(3)
m = vectorfy(M)
# h = (x[0]**2-x[0])*(x[1]**2-x[1])*(x[2]**2-x[2])*Matrix([1,1,1,1,1])

###

S0 = (const.B + sqrt(const.B**2 + 24.0*const.A*const.C))/(4.0*const.C)
Q0 = S0*(outerp(nu,nu) - (1.0/3.0)*eye(3))
QT = Q + S0*eye(3)
Pi = eye(3) - outerp(nu,nu)

# Energies

energyElastic = GeneralForm(const.L1/2*termL1(Dq,Dq)+const.L2/2*termL2(Dq,Dq)+const.L3/2*termL3(Dq,Dq),[Dq,q],name='energyElastic')
energyBulkC = GeneralForm((1/const.ep**2)*(((const.L0 - const.A)/2)*innerp(Q,Q) - (const.B/3)*trace(Q**3) + (const.C/4)*trace(Q**2)**2),[Dq,q])
energyBulkE = GeneralForm((1/const.ep**2)*(const.L0/2)*trace(Q**2),[Dq,q])

# Bilinear forms

bfTimeStep = GeneralForm((1/const.dt)*q.dot(p),[Dq,q],[Dp,p],name='bfTimeStep')
bfElastic = variationalDerivative(energyElastic,[Dq,q],[Dp,p],name='bfElastic')
bfBulkC = variationalDerivative(energyBulkC,[Dq,q],[Dp,p],name='bfBulkC')
bfBulkE = variationalDerivative(energyBulkE,[Dq,q],[Dp,p],name='bfBulkE')

bfTwist = GeneralForm(const.L1/2*term_twist_var(q,p),[Dq,q],[Dp,p],name='bfTwist')

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

class comp:
    initial_q = uflfy(m)
    manufac_q = uflfy(m)
    forcing_f = uflfy(vectorfy(strongForm(M)))
    forcing_g = uflfy(vectorfy(strongFormGamma(M)))

    n_bf_O = newt_bilinearDomain()
    n_bf_G = newt_bilinearBoundary()
    n_lf_O = newt_linearDomain()
    n_lf_G = newt_linearBoundary()

# END OF CODE