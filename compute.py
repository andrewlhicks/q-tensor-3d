from sympyplus import *
from settings import const, options, userfunc

# Set up variable with respect to which we will take derivatives

x = [Symbol('x0'),Symbol('x1'),Symbol('x2')]

n = Matrix([cos(x[0]+x[1]+x[2]),sin(x[0]+x[1]+x[2]),0])
G = outerp(n,n) - (1.0/3.0) * eye(3)
g = vectorfy(G)
# G = tensorfy(sin(N(2*pi)*x[0])*sin(N(2*pi)*x[1])*sin(N(2*pi)*x[2])*Matrix([1,1,1,1,1]))
# g = Matrix([(x[0]**2-x[0])*(x[1]**2-x[1])*(x[2]**2-x[2]),0,0,0,0])
# G = tensorfy(g)

# Intial guess and boundary

def manu_soln():
    return uflfy(g)

def manu_forc():
    return uflfy(vectorfy(strongForm(G)))

def bdy():
    return uflfy(g)

h = 0*(x[0]**2-x[0])*(x[1]**2-x[1])*(x[2]**2-x[2])*Matrix([1,1,1,1,1])
initial_guess = g + h

def initialGuess():
    # return uflfy(userfunc.initialGuess())
    return uflfy(initial_guess)

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
# energyNAnchor = GeneralForm(const.W0/2*innerp(Q - Q0,Q - Q0),[Dq,q])
# energyPDAnchor = GeneralForm(const.W1/2*innerp(QT-Pi*QT*Pi,QT-Pi*QT*Pi) + const.W2/2*(innerp(QT,QT) - S0**2)**2,[Dq,q])

# Bilinear forms

bfTimeStep = GeneralForm((1/const.dt)*q.dot(p),[Dq,q],[Dp,p],name='bfTimeStep')
bfElastic = variationalDerivative(energyElastic,[Dq,q],[Dp,p],name='bfElastic')
bfBulkC = variationalDerivative(energyBulkC,[Dq,q],[Dp,p],name='bfBulkC')
bfBulkE = variationalDerivative(energyBulkE,[Dq,q],[Dp,p],name='bfBulkE')
# bfPDAnchor = variationalDerivative(energyPDAnchor,[Dq,q],[Dp,p],name='bfPDAnchor')
# bfNAnchor = variationalDerivative(energyNAnchor,[Dq,q],[Dp,p],name='bfNAnchor')

# Linear forms

lfTimeStep = GeneralForm(bfTimeStep([Dqp,qp],[Dp,p]),[Dp,p],name='lfTimeStep')
lfBulkE = GeneralForm(bfBulkE([Dqp,qp],[Dp,p]),[Dp,p],name='lfBulkE')

# Lefthand side

bilinearDomain = lhsForm([Dq,q],[Dp,p],forms=[bfTimeStep,bfElastic,bfBulkC])

# Righthand side

linearDomain = rhsForm([Dp,p],forms=[lfTimeStep,lfBulkE])
if options.manufactured:
    lfManuForcing = GeneralForm(f.dot(p),[Dp,p],name='lfManuForcing')
    linearDomain.add_form(lfManuForcing)

# Newton's method

newt_bilinearDomain, newt_linearDomain = newtonsMethod(bilinearDomain,linearDomain,[Dqnp,qnp],[Dq,q],[Dp,p])

# END OF CODE