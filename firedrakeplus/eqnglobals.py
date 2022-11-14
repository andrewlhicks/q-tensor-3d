def set_eqn_globals(comp,uflcache_dict):
    from config import settings

    global EqnGlobals
    class EqnGlobals:
        """ It is evident that I shouldn't need to add to this class every time a variable is added to EqnGlobals.
        Something should be done about that. """
        # first, create three separate PDE systems
        pde = comp['pde']
        pde_nm = comp['pde_nm']
        pde_pd = comp['pde_pd']

        initial_q = uflcache_dict['initcond'] # uflcache will be preferred over comp in the future
        manu_q = uflcache_dict['manu_q'] if 'manu_q' in uflcache_dict.keys() else 'as_vector([0,0,0,0,0])'

        forcing_f = uflcache_dict['forcing_f'] if 'forcing_f' in uflcache_dict.keys() else 'as_vector([0,0,0,0,0])' # NOTE: while it works to calculate f here as an interpolation with Firedrake, it's actually more precise to do it in sympy beforehand.
        forcing_g = uflcache_dict['forcing_g'] if 'forcing_g' in uflcache_dict.keys() else 'as_vector([0,0,0,0,0])'
        s_bdy = uflcache_dict['s_bdy'] if 's_bdy' in uflcache_dict.keys() else 'as_vector([0,0,0,0,0])' # warning: may want to consider this carefully

        w_bdy_nu = uflcache_dict['w_bdy_nu'] if 'w_bdy_nu' in uflcache_dict.keys() else 'FacetNormal(mesh)'

        energies = comp['energies']
        energies_minmom = comp['energies_minmom']