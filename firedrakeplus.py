def errorH1(func_comp,func_true):
    from firedrake import sqrt, assemble, inner, grad, dot, dx
    diff = func_comp - func_true
    return sqrt(assemble((inner(grad(diff),grad(diff)) + dot(diff,diff)) * dx))

def errorL2(func_comp,func_true):
    from firedrake import sqrt, assemble, dot, dx
    diff = func_comp - func_true
    return sqrt(assemble((dot(diff,diff)) * dx))

def tensorfy(vector):
    from firedrake import sqrt, as_tensor
    
    # Basis of Q-tensor for Eigen calculation

    a = (sqrt(3.0)-3.0)/6.0
    b = (sqrt(3.0)+3.0)/6.0
    c = -sqrt(3.0)/3.0
    d = sqrt(2.0)/2.0

    E0 = as_tensor([[a,0,0],[0,b,0],[0,0,c]])
    E1 = as_tensor([[b,0,0],[0,a,0],[0,0,c]])
    E2 = as_tensor([[0,d,0],[d,0,0],[0,0,0]])
    E3 = as_tensor([[0,0,d],[0,0,0],[d,0,0]])
    E4 = as_tensor([[0,0,0],[0,0,d],[0,d,0]])
    
    # Return the linear combination of tensors
    
    return vector[0] * E0 + vector[1] * E1 + vector[2] * E2 + vector[3] * E3 + vector[4] * E4

# END OF CODE