def valueCheck():
    from settings import const, options
    from misc import color

    L1, L2, L3 = const['L1'], const['L2'], const['L3']
    
    if not isinstance(options['manufactured'],bool):
        raise ValueError("Variable 'manufactured' must be a boolean.")

    if not isinstance(options['omit_init_printoff'],bool):
        raise ValueError("Variable 'omit_init_printoff' must be a boolean.")

    if not isinstance(options['visualize'],bool):
        raise ValueError("Variable 'visualize' must be a boolean.")
    
    if 0>=L1 or -L1>=L3 or L3>=2*L1 or -3/5*L1-1/10*L3>=L2:
        print()
        print(f"{color.warning}WARNING: L1, L2, and L3 do not satisfy the proper inequalities{color.end}")