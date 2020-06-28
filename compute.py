import settings

def computeBilinear():
    from settings import dt, ep, L0, L1, L2, L3
    from sympy import symbols, Matrix
    from sympyplus import uflfy
    
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
    
    qv = Matrix([q0,q1,q2,q3,q4])
    pv = Matrix([p0,p1,p2,p3,p4])
    
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
    
    # Computations
    
    term_L1 = termL1(grad_q,grad_p)
    term_L2 = termL2(grad_q,grad_p)
    term_L3 = termL3(grad_q,grad_p)
    
    # Combine
    
    expression = (1/dt) * qv.dot(pv) + L1*term_L1 + L2*term_L2 + L3*term_L3 + (L0/ep)*qv.dot(pv)
    
    return uflfy(expression)

def computeBoundary():
    from sympy import eye
    from sympyplus import outerp, uflfy, vectorfy
    
    n = settings.boundary()
    G = outerp(n,n) - (1.0/3.0) * eye(3)
    
    return uflfy(vectorfy(G))

def computeInitialGuess():
    from sympyplus import uflfy
    return uflfy(settings.initialGuess())

def computeLinear():
    from settings import dt, ep, L0, A, B, C
    from sympy import symbols, sqrt, diff, diag, Matrix, zeros, eye
    from sympyplus import uflfy, innerp, outerp, vectorfy, E
    
    qp0, qp1, qp2, qp3, qp4 = symbols('q_prev[0:5]')
    qpv = Matrix([qp0,qp1,qp2,qp3,qp4])
    
    p0, p1, p2, p3, p4 = symbols('p[0:5]')
    pv = Matrix([p0,p1,p2,p3,p4])
    
    fv = strongForm()
        
    term_A = termA(qpv,pv)
    term_B = termB(qpv,pv)
    term_C = termC(qpv,pv)
    
    return uflfy((1/dt) * qpv.dot(pv) + (1/ep) * ( (A + L0) * term_A + B * term_B - C * term_C ) + fv.dot(pv))

######################################################################################################################################################################################
######################################################################################################################################################################################
######################################################################################################################################################################################

def strongForm():
    from settings import ep, L1, L2, L3, A, B, C
    from sympy import Matrix, eye
    from sympyplus import outerp, vectorfy
    
    if settings.manufactured == 0:
        return Matrix([0,0,0,0,0])
    elif settings.manufactured == 1:
        n = settings.boundary()
        G = outerp(n,n) - (1.0/3.0) * eye(3)
        
        strong_L1 = strongL1(G)
        strong_L2 = strongL2(G)
        strong_L3 = strongL3(G)
        strong_A = strongA(G)
        strong_B = strongB(G)
        strong_C = strongC(G)
        
        F = L1*strong_L1 + L2*strong_L2 + L3*strong_L3 + (1/ep)*(-A*strong_A - B*strong_B + C*strong_C)
        
        return vectorfy(F)
    else:
        raise ValueError("Variable 'manufactured' must be 0 or 1.")

def termA(vec1,vec2):
    return vec1.dot(vec2)

def termB(vec1,vec2):
    from sympy import zeros
    from sympyplus import E, innerp
    
    ten1 = zeros(3,3)
    ten2 = zeros(3,3)
    
    for ii in range(5):
        ten1 += vec1[ii]*E[ii]
        ten2 += vec2[ii]*E[ii]
    
    return innerp(ten1*ten1,ten2)

def termC(vec1,vec2):
    return vec1.dot(vec1) * vec1.dot(vec2)

def termL1(grad1,grad2):
    from sympyplus import innerp
    
    return innerp(grad1,grad2)

def termL2(grad1,grad2):
    from sympyplus import E
    
    term = 0

    for ii in range(5):
        for jj in range(5):
            term += (grad1.row(ii) * E[ii]).dot(E[jj] * grad2.row(jj).T)
    
    return term

def termL3(grad1,grad2):
    from sympyplus import E
    
    term = 0

    for ii in range(5):
        for jj in range(5):
            term += (grad1.row(ii) * E[jj]).dot(E[ii] * grad2.row(jj).T)
    
    return term

######################################################################################################################################################################################
######################################################################################################################################################################################
######################################################################################################################################################################################

def strongL1(Q):
    from sympy import symbols, zeros, diff
    x0,x1,x2 = symbols('x0 x1 x2')
    x = [x0,x1,x2]
    
    term = zeros(3,3)
    
    for ii in range(3):
            for jj in range(3):
                for kk in range(3):
                    term[ii,jj] -= diff(Q[ii,jj],x[kk],x[kk])
    
    return term

def strongL2(Q):
    from sympy import symbols, zeros, diff
    x0,x1,x2 = symbols('x0 x1 x2')
    x = [x0,x1,x2]
    
    term = zeros(3,3)
        
    for ii in range(3):
        for jj in range(3):
            for kk in range(3):
                term[ii,jj] -= diff(Q[ii,kk],x[jj],x[kk])
    
    return term

def strongL3(Q):
    from sympy import symbols, zeros, diff
    x0,x1,x2 = symbols('x0 x1 x2')
    x = [x0,x1,x2]
    
    term = zeros(3,3)
        
    for ii in range(3):
        for jj in range(3):
            for kk in range(3):
                term[ii,jj] -= diff(Q[ii,kk],x[kk],x[jj])
    
    return term

def strongA(Q):
    return Q

def strongB(Q):
    return Q*Q

def strongC(Q):
    from sympyplus import innerp
    
    return innerp(Q,Q)*Q

# END OF CODE