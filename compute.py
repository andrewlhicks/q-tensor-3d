from sympyplus import *
from misc import getValues
from settings import const, options, userfunc

# Set up variable with respect to which we will take derivatives

x = [Symbol('x0'),Symbol('x1'),Symbol('x2')]

# Intial guess and boundary

def boundary():   
    n = userfunc.boundary()
    G = outerp(n,n) - (1.0/3.0) * eye(3)
    
    return uflfy(vectorfy(G))

def initialGuess():
    return uflfy(userfunc.initialGuess())

# Treatment of the L1, L2, and L3 terms, since they are too complicated to write in-line

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

def tensorfyFromUnitVec(unit_vec):
    return outerp(unit_vec,unit_vec) - (1.0/3.0) * eye(3)

def strongForm(G): # plugs G into the strong form PDE
    return const.L1*strongL1(G) + const.L2*strongL2(G) + const.L3*strongL3(G) + (1/const.ep**2)*(-const.A*G - const.B*G*G + const.C*innerp(G,G)*G)

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

qnp = QVector('q_newt_prev')
Dqnp = qnp.grad
QNP = qnp.tens

##############

energyElastic = GeneralForm(const.L1/2*termL1(Dq,Dq)+const.L2/2*termL2(Dq,Dq)+const.L3/2*termL3(Dq,Dq),[Dq,q],name='energyElastic')

# New convex splitting
# energyBulkC = GeneralForm((1/const.ep**2)*((const.L0/2 - const.A)*innerp(Q,Q) - (const.B/3)*trace(Q**3) + (const.C/4)*trace(Q**2)**2),[Dq,q])
# energyBulkE = GeneralForm((1/const.ep**2)*(const.L0/2)*trace(Q**2),[Dq,q])

# Old convex splitting
energyBulkC = GeneralForm((1/const.ep**2)*(const.L0/2)*trace(Q**2),[Dq,q],name='energyBulkC')
energyBulkE = GeneralForm((1/const.ep**2)*((const.L0+const.A)/2*trace(Q**2) + (const.B/3)*trace(Q**3) - (const.C/4)*trace(Q**2)**2),[Dq,q],name='energyBulkE')

W0 = 1
W1 = 1
W2 = 1

S0 = (const.B + sqrt(const.B**2 + 24.0*const.A*const.C))/(4.0*const.C)
Q0 = S0*(outerp(nu,nu) - (1.0/3.0)*eye(3))
QT = Q + S0*eye(3)
Pi = eye(3) - outerp(nu,nu)

energyNAnchor = GeneralForm(W0*innerp(Q - Q0,Q - Q0),[Dq,q])
energyPDAnchor = GeneralForm(W1*innerp(QT-Pi*QT*Pi,QT-Pi*QT*Pi) + W2*(innerp(QT,QT) - S0**2)**2,[Dq,q])

##############

bfTimeStep = GeneralForm((1/const.dt)*q.dot(p),[Dq,q],[Dp,p],name='bfTimeStep')
bfElastic = variationalDerivative(energyElastic,[Dq,q],[Dp,p],name='bfElastic')
bfBulkC = variationalDerivative(energyBulkC,[Dq,q],[Dp,p],name='bfBulkC')
bfBulkE = variationalDerivative(energyBulkE,[Dq,q],[Dp,p],name='bfBulkE')
bfPDAnchor = variationalDerivative(energyPDAnchor,[Dq,q],[Dp,p],name='bfPDAnchor')
bfNAnchor = variationalDerivative(energyNAnchor,[Dq,q],[Dp,p],name='bfNAnchor')

##############

lfTimeStep = GeneralForm(bfTimeStep([Dqp,qp],[Dp,p]),[Dp,p],name='lfTimeStep')
lfBulkE = GeneralForm(bfBulkE([Dqp,qp],[Dp,p]),[Dp,p],name='lfBulkE')
lfManuForcing = GeneralForm(vectorfy(strongForm(tensorfyFromUnitVec(userfunc.boundary()))).dot(p),[Dp,p],name='lfManuForcing')

##############

# bilinearDomain = BilinearForm([Dq,q],[Dp,p])
# bilinearDomain.add_form((1/const.dt) * q.dot(p))
# bilinearDomain.add_form(tfElastic(Dqnp,qnp,[Dq,q],[Dp,p]))
# bilinearDomain.add_form((1/const.ep**2) * tfBulkC[(Dqnp,qnp],[Dq,q],[Dp,p]))

# linearDomain = LinearForm([Dp,p)]
# # First subtract the original a(q,p)
# linearDomain.add_form(- (1/const.dt) * qnp.dot(p))
# linearDomain.add_form(- bfElastic([Dqnp,qnp],[Dp,p]))
# linearDomain.add_form(- (1/const.ep**2) * bfBulkC([Dqnp,qnp],[Dp,p]))
# # Then add the original L(p)
# linearDomain.add_form((1/const.ep**2) * bfBulkE([Dqp,qp],[Dp,p]))
# linearDomain.add_form((1/const.dt) * qp.dot(p))
# linearDomain.add_form(lfForcing([Dp,p]))

bilinearDomain = lhsForm([Dq,q],[Dp,p],forms=[bfTimeStep,bfElastic,bfBulkC])

linearDomain = rhsForm([Dp,p],forms=[lfTimeStep,lfBulkE])
if options.manufactured: linearDomain.add_form(lfManuForcing)

# bilinearBoundary = BilinearForm([Dq,q],[Dp,p])
# bilinearBoundary.add_form(tfNAnchor([Dqnp,qnp],[Dq,q],[Dp,p]))
# bilinearBoundary.add_form(tfPDAnchor([Dqnp,qnp],[Dq,q],[Dp,p]))
# linearBoundary = LinearForm([Dp,p])
# linearBoundary.add_form(- bfNAnchor([Dqnp,qnp],[Dp,p]))
# linearBoundary.add_form(- bfPDAnchor([Dqnp,qnp],[Dp,p]))

newt_bilinearDomain, newt_linearDomain = newtonsMethod(bilinearDomain,linearDomain,[Dqnp,qnp],[Dq,q],[Dp,p])

# END OF CODE