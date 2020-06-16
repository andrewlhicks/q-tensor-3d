from sympy import *

def inner(A,B):
    return trace(A.T*B)

def outer(A,B):
    return A*B.T

a = (sqrt(3)-3)/6
b = (sqrt(3)+3)/6
c = -sqrt(3)/3
d = sqrt(2)/2

E = [diag(a,b,c), diag(b,a,c), Matrix([[0,d,0],[d,0,0],[0,0,0]]), Matrix([[0,0,d],[0,0,0],[d,0,0]]), Matrix([[0,0,0],[0,0,d],[0,d,0]])]


L0,ep,dt,L1,L2,L3,A,B,C = symbols('L0,ep,dt,L1,L2,L3,A,B,C')

# L0 = 10.0
# ep = 1000000
# dt = 0.1

# L1 = 1
# L2 = 0
# L3 = 0

# A = 2
# B = 0
# C = 4

x,y,z = symbols('x0 x1 x2')

n = Matrix([sin(x),cos(y),1])
G = outer(n,n) - (1/3.0) * eye(3)
gv = Matrix([inner(G,E[0]),inner(G,E[1]),inner(G,E[2]),inner(G,E[3]),inner(G,E[4])])

grad_g = Matrix([[diff(gv[0],x),diff(gv[0],y),diff(gv[0],z)],
                     [diff(gv[1],x),diff(gv[1],y),diff(gv[1],z)],
                     [diff(gv[2],x),diff(gv[2],y),diff(gv[2],z)],
                     [diff(gv[3],x),diff(gv[3],y),diff(gv[3],z)],
                     [diff(gv[4],x),diff(gv[4],y),diff(gv[4],z)]])

strong_L1 = diff(grad_g,x) + diff(grad_g,y) + diff(grad_g,z)

Mat = Matrix([[0,1],[-1,0]])
# fv = -(L1+(L2+L3)/2)*(diff(gv,x,2)+diff(gv,y,2)) + (L2-L3)/2*Mat*(diff(gv,x,y)-diff(gv,y,x)) - (2/ep)*gv + (4/ep)*gv.dot(gv)*gv
fv = Matrix([0,0,0,0,0])

q0, q1, q2, q3, q4 = symbols('q[0:5]')

d0q0, d1q0, d2q0, \
d0q1, d1q1, d2q1, \
d0q2, d1q2, d2q2, \
d0q3, d1q3, d2q3, \
d0q4, d1q4, d2q4 = symbols('q[0:5].dx((0:3))')

p0, p1, p2, p3, p4 = symbols('p[0:5]')

d0p0, d1p0, d2p0, \
d0p1, d1p1, d2p1, \
d0p2, d1p2, d2p2, \
d0p3, d1p3, d2p3, \
d0p4, d1p4, d2p4 = symbols('p[0:5].dx((0:3))')


qp0, qp1, qp2, qp3, qp4 = symbols('q_prev[0:5]')

qv = Matrix([q0,q1,q2,q3,q4])
pv = Matrix([p0,p1,p2,p3,p4])

qpv = Matrix([qp0,qp1,qp2,qp3,qp4])

grad_q = Matrix([[d0q0,d1q0,d2q0],
                     [d0q1,d1q1,d2q1],
                     [d0q2,d1q2,d2q2],
                     [d0q3,d1q3,d2q3],
                     [d0q4,d1q4,d2q4]])
                     
grad_p = Matrix([[d0p0,d1p0,d2p0],
                     [d0p1,d1p1,d2p1],
                     [d0p2,d1p2,d2p2],
                     [d0p3,d1p3,d2p3],
                     [d0p4,d1p4,d2p4]])

term_L1 = inner(grad_q,grad_p)

term_L2 = 0

for ii in range(5):
    for jj in range(5):
        term_L2 += (grad_q.row(ii) * E[ii]).dot(E[jj] * grad_p.row(jj).T)

term_L3 = 0

for ii in range(5):
    for jj in range(5):
        term_L3 += (grad_q.row(ii) * E[jj]).dot(E[ii] * grad_p.row(jj).T)

term_A = qpv.dot(pv)

# QP:Q + PQ:Q + Q^2:P


Qt = zeros(3,3)
Pt = zeros(3,3)

for ii in range(5):
    Qt += qv[ii]*E[ii]
    Pt += pv[ii]*E[ii]

Qt = zeros(2,2)
Pt = zeros(2,2)

F = [1/sqrt(2)*Matrix([[1,0],[0,-1]]), 1/sqrt(2)*Matrix([[0,1],[1,0]])]

for ii in range(2):
    Qt += qv[ii]*F[ii]
    Pt += pv[ii]*F[ii]

term_B = (inner(Qt*Pt,Qt) + inner(Pt*Qt,Qt) + inner(Qt*Qt,Pt) )/3.0

term_C = qpv.dot(qpv) * qpv.dot(pv)

a0,a1,a2,a3 = symbols('a0:4')
b0,b1,b2,b3 = symbols('b0:4')
c0,c1,c2,c3 = symbols('c0:4')

A = Matrix([[a0,a1],[a2,a3]])
B = Matrix([[b0,b1],[b2,b3]])
C = Matrix([[c0,c1],[c2,c3]])

# bilinear_a = (1/dt) * qv.dot(pv) + L1*term_L1 + L2*term_L2 + L3*term_L3 + (L0/ep)*qv.dot(pv)
# linear_L   = (1/dt) * qpv.dot(pv) + (1/ep) * ( (A + L0) * term_A - B * term_B - C * term_C ) + fv.dot(pv)