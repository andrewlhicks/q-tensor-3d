def set_eqn_globals(comp,uflcache_dict):
    # NOTES: (1) While it works to calculate f here as an interpolation with Firedrake, it's actually more precise to do it in sympy beforehand.
    #        (2) We may want to consider carefully how we handle the default state of 's_bdy'
    #        (3) It is evident that I shouldn't need to add to this class every time a variable is added to EqnGlobals. Something should be done about that.
    global EqnGlobals
    class EqnGlobals:
        # comp: PDE systems and energies
        pde = comp['pde']
        pde_nm = comp['pde_nm']
        pde_pd = comp['pde_pd']
        energies = comp['energies']
        energies_minmom = comp['energies_minmom']

        # uflcache dicts: will be preferred over comp in the future
        process_uflcache_dict(uflcache_dict)
        initial_q = uflcache_dict['initcond'] 
        manu_q = uflcache_dict['manu_q']
        forcing_f = uflcache_dict['forcing_f']
        forcing_g = uflcache_dict['forcing_g']
        s_bdy = uflcache_dict['s_bdy']
        w_bdy_nu = uflcache_dict['w_bdy_nu']

def process_uflcache_dict(uflcache_dict):
    # if the 'initcond' key is missing, raise an error
    if 'initcond' not in uflcache_dict.keys(): raise ValueError('Missing "initcond" in ulfcache')
    # if the following four keys are missing, then add the key and make its value the zero vector in 5d
    [uflcache_dict.update({key:'as_vector([0,0,0,0,0])'}) for key in ('manu_q','forcing_f','forcing_g','s_bdy') if key not in uflcache_dict.keys()]
    # if the 'w_bdy_nu' key is missing, then add the key and make its value the facet normal from the mesh
    if 'w_bdy_nu' not in uflcache_dict.keys(): uflcache_dict.update({'w_bdy_nu':'FacetNormal(mesh)'})