def set_eqn_globals(comp,uflcache_dict):
    from config import settings

    global EqnGlobals
    class EqnGlobals:
        # first, create three separate PDE systems
        pde = comp['pde']
        pde_nm = comp['pde_nm']
        pde_pd = comp['pde_pd']

        initial_q = uflcache_dict['initcond'] # uflcache will be preferred over comp in the future

        forcing_f = comp['forcing_f'] if settings.options.manufactured else 'as_vector([0,0,0,0,0])' # NOTE: while it works to calculate f here as an interpolation with Firedrake, it's actually more precise to do it in sympy beforehand.
        forcing_g = comp['forcing_g'] if settings.options.manufactured else 'as_vector([0,0,0,0,0])'
        bdy_cond = comp['bdy_cond']

        energies = comp['energies']