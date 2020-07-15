def valueCheck():
    from settings import L1, L2, L3, manufactured, omit_init_printoff, visualize
    from misc import color

    if manufactured != 0 and manufactured != 1:
        raise ValueError("Variable 'manufactured' must be set to 0 or 1.")

    if omit_init_printoff != 0 and omit_init_printoff != 1:
        raise ValueError("Variable 'omit_init_printoff' must be set to 0 or 1.")

    if visualize != 0 and visualize != 1:
        raise ValueError("Variable 'visualize' must be set to 0 or 1.")
    
    if 0>=L1 or -L1>=L3 or L3>=2*L1 or -3/5*L1-1/10*L3>=L2:
        print()
        print(f"{color.warning}WARNING: L1, L2, and L3 do not satisfy the proper inequalities{color.end}")