from sympyplus import *
from compute_terms import *
from settings import const, options

# Initial guess and manufactured solution

n = Matrix([cos(x[0]+x[1]+x[2]),sin(x[0]+x[1]+x[2]),0])
G = outerp(n,n) - (1.0/3.0) * eye(3)
g = vectorfy(G)
h = 0*(x[0]**2-x[0])*(x[1]**2-x[1])*(x[2]**2-x[2])*Matrix([1,1,1,1,1])

# Intial guess and boundary

def init_guess():
    return uflfy(g + h)

def manu_soln():
    return uflfy(g)

def manu_forc():
    return uflfy(vectorfy(strongForm(G)))

# Set up Qvector objects

nu = AbstractVector('nu')

q = QVector('q')
Dq = q.grad
Q = q.tens

p = QVector('p')
Dp = p.grad
P = p.tens

qp = QVector('q_prev')
Dqp = qp.grad
QP = qp.tens

qnp = QVector('q_newt_prev')
Dqnp = qnp.grad
QNP = qnp.tens

f = QVector('f')

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

bfNAnchor = GeneralForm(const.W0*q.dot(p),[Dq,q],[Dp,p],name='bfNAnchor')
bfPDAnchor1 = GeneralForm(const.W1*innerp(Q-Pi*Q*Pi,P-Pi*P*Pi),[Dq,q],[Dp,p],name='bfPDAnchor1')
bfPDAnchor2 = GeneralForm(const.W2*(innerp(Q,Q)*innerp(Q,P)-8*S0**2/9*innerp(Q,P)),[Dq,q],[Dp,p],name='bfPDAnchor2')

# Linear forms

lfTimeStep = GeneralForm(bfTimeStep([Dqp,qp],[Dp,p]),[Dp,p],name='lfTimeStep')
lfBulkE = GeneralForm(bfBulkE([Dqp,qp],[Dp,p]),[Dp,p],name='lfBulkE')

lfNAnchor = GeneralForm(const.W0*innerp(Q0,P),[Dp,p],name='lfNAnchor')
lfPDAnchor1 = GeneralForm(const.W1*S0/3.0*innerp(eye(3)-Pi,P-Pi*P*Pi),[Dp,p],name='lfPDAnchor1')

# Lefthand side

bilinearDomain = lhsForm([Dq,q],[Dp,p],forms=[bfTimeStep,bfElastic,bfBulkC])
bilinearBoundary = lhsForm([Dq,q],[Dp,p],forms=[bfNAnchor,bfPDAnchor1,bfPDAnchor2])

# Righthand side

linearDomain = rhsForm([Dp,p],forms=[lfTimeStep,lfBulkE])
if options.manufactured:
    lfManuForcing = GeneralForm(f.dot(p),[Dp,p],name='lfManuForcing')
    linearDomain.add_form(lfManuForcing)
linearBoundary = rhsForm([Dp,p],forms=[lfNAnchor,lfPDAnchor1])

# Newton's method

newt_bilinearDomain, newt_linearDomain = newtonsMethod(bilinearDomain,linearDomain,[Dqnp,qnp],[Dq,q],[Dp,p])
newt_bilinearBoundary, newt_linearBoundary = newtonsMethod(bilinearBoundary,linearBoundary,[Dqnp,qnp],[Dq,q],[Dp,p])

# Create relevant UFL strings

class comp:
    init_guess = init_guess()
    manu_soln = manu_soln()
    manu_forc = manu_forc()

    n_bf_O = newt_bilinearDomain()
    n_bf_G = newt_bilinearBoundary()
    n_lf_O = newt_linearDomain()
    n_lf_G = newt_linearBoundary()

# END OF CODE