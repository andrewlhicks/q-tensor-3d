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

r = QVector('r')
Dr = r.grad
R = r.tens

qp = QVector('q_prev')
Dqp = qp.grad
QP = qp.tens

qnp = QVector('q_newt_prev')
Dqnp = qnp.grad
QNP = qnp.tens

f = QVector('f')

##############

energyElastic = GeneralForm(const.L1/2*termL1(Dq,Dq)+const.L2/2*termL2(Dq,Dq)+const.L3/2*termL3(Dq,Dq),[Dq,q],name='energyElastic')

# New convex splitting
energyBulkC = GeneralForm((1/const.ep**2)*(((const.L0 - const.A)/2)*innerp(Q,Q) - (const.B/3)*trace(Q**3) + (const.C/4)*trace(Q**2)**2),[Dq,q])
energyBulkE = GeneralForm((1/const.ep**2)*(const.L0/2)*trace(Q**2),[Dq,q])

# Old convex splitting
# energyBulkC = GeneralForm((1/const.ep**2)*(const.L0/2)*trace(Q**2),[Dq,q],name='energyBulkC')
# energyBulkE = GeneralForm((1/const.ep**2)*((const.L0+const.A)/2*trace(Q**2) + (const.B/3)*trace(Q**3) - (const.C/4)*trace(Q**2)**2),[Dq,q],name='energyBulkE')

# W0 = 1
# W1 = 1
# W2 = 1

# S0 = (const.B + sqrt(const.B**2 + 24.0*const.A*const.C))/(4.0*const.C)
# Q0 = S0*(outerp(nu,nu) - (1.0/3.0)*eye(3))
# QT = Q + S0*eye(3)
# Pi = eye(3) - outerp(nu,nu)

# energyNAnchor = GeneralForm(W0*innerp(Q - Q0,Q - Q0),[Dq,q])
# energyPDAnchor = GeneralForm(W1*innerp(QT-Pi*QT*Pi,QT-Pi*QT*Pi) + W2*(innerp(QT,QT) - S0**2)**2,[Dq,q])

##############

bfTimeStep = GeneralForm((1/const.dt)*q.dot(p),[Dq,q],[Dp,p],name='bfTimeStep')
bfElastic = variationalDerivative(energyElastic,[Dq,q],[Dp,p],name='bfElastic')
bfBulkC = variationalDerivative(energyBulkC,[Dq,q],[Dp,p],name='bfBulkC')
bfBulkE = variationalDerivative(energyBulkE,[Dq,q],[Dp,p],name='bfBulkE')
# bfPDAnchor = variationalDerivative(energyPDAnchor,[Dq,q],[Dp,p],name='bfPDAnchor')
# bfNAnchor = variationalDerivative(energyNAnchor,[Dq,q],[Dp,p],name='bfNAnchor')

##############

lfTimeStep = GeneralForm(bfTimeStep([Dqp,qp],[Dp,p]),[Dp,p],name='lfTimeStep')
lfBulkE = GeneralForm(bfBulkE([Dqp,qp],[Dp,p]),[Dp,p],name='lfBulkE')


lfManuForcing = GeneralForm(f.dot(p),[Dp,p],name='lfManuForcing')

##############

# der_bfBulkC = GeneralForm( (1/const.ep**2) * ( (const.L0 - const.A)*innerp(R,P) - 2*const.B*innerp(Q*R,P) + const.C*(2*innerp(Q,R)*innerp(Q,P) + innerp(Q,Q)*innerp(R,P) ) ), [Dq,q],[Dr,r],[Dp,p])

##############

bilinearDomain = lhsForm([Dq,q],[Dp,p],forms=[bfTimeStep,bfElastic,bfBulkC])

linearDomain = rhsForm([Dp,p],forms=[lfTimeStep,lfBulkE])
if options.manufactured: linearDomain.add_form(lfManuForcing)

newt_bilinearDomain, newt_linearDomain = newtonsMethod(bilinearDomain,linearDomain,[Dqnp,qnp],[Dq,q],[Dp,p])

##############

# newt_bilinearDomain = lhsForm([Dq,q],[Dp,p])
# newt_bilinearDomain.add_form(bfTimeStep.new_params([Dq,q],[Dp,p]))
# newt_bilinearDomain.add_form(bfElastic.new_params([Dq,q],[Dp,p]))
# newt_bilinearDomain.add_form(der_bfBulkC.new_params([Dqnp,qnp],[Dq,q],[Dp,p]))

# newt_linearDomain = rhsForm([Dp,p])
# # First subtract the original a(q,p)
# newt_linearDomain.add_form(bfTimeStep.new_params([Dqnp,qnp],[Dp,p]).mul(-1))
# newt_linearDomain.add_form(bfElastic.new_params([Dqnp,qnp],[Dp,p]).mul(-1))
# newt_linearDomain.add_form(bfBulkC.new_params([Dqnp,qnp],[Dp,p]).mul(-1))

# # Then add the original L(p)
# newt_linearDomain.add_form(lfBulkE)
# newt_linearDomain.add_form(lfTimeStep)
# if options.manufactured: newt_linearDomain.add_form(lfManuForcing)

# END OF CODE