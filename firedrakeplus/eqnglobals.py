from config import settings

def set_eqn_globals(comp,uflcache_dict):
    global EqnGlobals
    class EqnGlobals:
        # remove these next commit
        bilinear_form = comp['pde_d']['lhs']
        bilinear_form_bdy = comp['pde_b']['lhs']
        linear_form = comp['pde_d']['rhs']
        linear_form_bdy = comp['pde_b']['rhs']

        # make these standard next commit
        pde_d = comp['pde_d']
        pde_b = comp['pde_b']

        initial_q = uflcache_dict['initcond'] # uflcache will be preferred over comp in the future

        forcing_f = comp['forcing_f'] if settings.options.manufactured else 'as_vector([0,0,0,0,0])' # NOTE: while it works to calculate f here as an interpolation with Firedrake, it's actually more precise to do it in sympy beforehand.
        forcing_g = comp['forcing_g'] if settings.options.manufactured else 'as_vector([0,0,0,0,0])'
        strong_boundary = [comp['bdycond_s'],settings.options.strong_boundary]
        weak_boundary = [None,settings.options.weak_boundary]

        energies = comp['energies']