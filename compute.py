import settings

# Compute boundary and initial guess

def boundary():
    from settings import userBoundary
    from sympy import eye
    from sympyplus import outerp, uflfy, vectorfy
    
    n = userBoundary()
    G = outerp(n,n) - (1.0/3.0) * eye(3)
    
    return uflfy(vectorfy(G))

def initialGuess():
    from settings import userInitialGuess
    from sympyplus import uflfy
    
    return uflfy(userInitialGuess())

#  Compute bilinear forms integrated over the volume and on the boundary, then the linear form

def bilinear():
    from sympyplus import uflfy
    
    # Combine
        
    expression = bFormTimeStep() + bFormElastic() + bFormBulk()
    
    # Convert to UFL and return
    
    return uflfy(expression)

def bilinearOnBoundary():
    from sympyplus import uflfy
    
    # Combine
    
    expression = bFormSurface()
    
    # Convert to UFL and return
    
    return uflfy(expression)

def linear():
    from sympyplus import uflfy
    
    # Combine
    
    expression = lFormTimeStep() + lFormBulk() + lFormForcing()
    
    # Convert to UFL and return
    
    return uflfy(expression)

def linearOnBoundary():
    from sympyplus import uflfy
    
    # Combine
    
    expression = lFormSurface()
    
    # Convert to UFL and return
    
    return uflfy(expression)

######################################################################################################################################################################################
######################################################################################################################################################################################
######################################################################################################################################################################################

# BILINEAR FORMS

def bFormElastic():
    from settings import L1, L2, L3
    from sympyplus import uflfy
    
    # Create vector objects
    
    q = Vector('q')
    p = Vector('p')
    
    # Combine and return
    
    return L1*termL1(q.grad,p.grad) + L2*termL2(q.grad,p.grad) + L3*termL3(q.grad,p.grad)

def bFormBulk():
    from settings import dt, ep, L0
    from sympyplus import uflfy
    
    # Create vector objects
    
    q = Vector('q')
    p = Vector('p')
    
    # Combine and return
        
    return (L0/ep)*q.vec.dot(p.vec)

def bFormSurface():    
    # Create vector objects
    
    q = Vector('q')
    p = Vector('p')
    
    # Combine and return
    
    return q.vec.dot(p.vec)

def bFormTimeStep():
    from settings import dt
    
    # Create vector objects
    
    q = Vector('q')
    p = Vector('p')
    
    # Combine and return
    
    return (1/dt) * q.vec.dot(p.vec)

# LINEAR FORMS

def lFormBulk():
    from settings import ep, L0, A, B, C
    from sympyplus import uflfy
    
    # Create vector objects
    
    q_prev = Vector('q_prev')
    p = Vector('p')
    
    # Combine and return
    
    return (1/ep) * ( (A + L0) * termA(q_prev.vec,p.vec) + B * termB(q_prev.vec,p.vec) - C * termC(q_prev.vec,p.vec) )

def lFormForcing():
    from settings import userBoundary, manufactured
    from sympy import Matrix, eye
    from sympyplus import uflfy, outerp, vectorfy
    
    # Create vector object
    
    p = Vector('p')
    
    # Create a forcing if the manufactured solution is set
    
    if manufactured:
        n = userBoundary()
        G = outerp(n,n) - (1.0/3.0) * eye(3)
        F = strongForm(G)
        f = vectorfy(F)
    else:
        f = Matrix([0,0,0,0,0])
    
    # Combine and return
    
    return f.dot(p.vec)

def lFormSurface():
    from settings import userBoundary
    from sympy import eye
    from sympyplus import vectorfy, outerp
    
    # Create vector object
    
    p = Vector('p')
    
    # Compute Q_0 boundary tensor and q_0 vector
    
    # n = userBoundary()
    # Q_0 = outerp(n,n) - (1.0/3.0) * eye(3)
    # q_0 = vectorfy(Q_0)
    
    nu = Vector3d('nu')
    Q_0 = outerp(nu.vec,nu.vec) - (1.0/3.0) * eye(3) # Should be multiplied by s_0 value (min point for double well)
    q_0 = vectorfy(Q_0)

    # Combine and return
    
    return q_0.dot(p.vec)

def lFormTimeStep():
    from settings import dt
    from sympyplus import uflfy
    
    # Create vector object
    
    q_prev = Vector('q_prev')
    p = Vector('p')
    
    # Combine and return
    
    return (1/dt) * q_prev.vec.dot(p.vec)

######################################################################################################################################################################################
######################################################################################################################################################################################
######################################################################################################################################################################################

def strongForm(G): # plugs G into the strong form PDE
    from settings import ep, L1, L2, L3, A, B, C
    
    return L1*strongL1(G) + L2*strongL2(G) + L3*strongL3(G) + (1/ep)*(-A*strongA(G) - B*strongB(G) + C*strongC(G))

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

######################################################################################################################################################################################
######################################################################################################################################################################################
######################################################################################################################################################################################

class Vector: # creates a Vector object, which symbolically represents a 5-dimensional vector in 3 spatial dimensions, and its gradient matrix
    def __init__(self,name):
        self.name = name # For 'name', choose the variable name that Firedrake will later use
        
        from sympy import symbols, Matrix        
        v0, v1, v2, v3, v4 = symbols(f'{name}[0:5]')
        
        d0v0, d1v0, d2v0, \
        d0v1, d1v1, d2v1, \
        d0v2, d1v2, d2v2, \
        d0v3, d1v3, d2v3, \
        d0v4, d1v4, d2v4 = symbols(f'{name}[0:5].dx((0:3))')
        
        self.vec = Matrix([v0,v1,v2,v3,v4]) # the vector itself
        
        self.grad = Matrix([[d0v0,d1v0,d2v0],
                            [d0v1,d1v1,d2v1],
                            [d0v2,d1v2,d2v2],
                            [d0v3,d1v3,d2v3],
                            [d0v4,d1v4,d2v4]]) # the vector's gradient matrix

class Vector3d: # Temporary class until I can generalize the Vector class to handle this case
    def __init__(self,name):
        self.name = name
        from sympy import symbols, Matrix
        v0, v1, v2 = symbols(f'{name}[0:3]')
        self.vec = Matrix([v0,v1,v2])

# END OF CODE