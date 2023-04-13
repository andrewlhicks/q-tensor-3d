# Help

## Saves

To use qtensor3d, a save folder needs to be created. The easiest way to do this is to use `qsave` in the following manner:
```
qsave -b <savepath>
```
where `<savepath>` is the save folder. This will create a save folder with the default settings, cosntants, and uflexpr files.

## Settings configuration

All settings files are contained in the folder 'settings'. Settings files are in the '.yml' format. Visit [PyYAML](https://pyyaml.org/wiki/PyYAMLDocumentation) for more information.

- `mesh`:
  - `source` - `builtin`, `local`, `global`, or `legacy`
  - `name` - See "Choosing a mesh" for more info
  - `radius` - The characteristic length for non-dimensionalization
  - `refs` - How many refinements of the mesh to complete
- `options`:
  - `nondim` - If true, nondimensionalizes the constants before simulation
  - `strong_boundary` - May be given as a number or a list of numbers; specifies those boundaries where Dirichlet conditions will be enforced
  - `weak_boundary` - May be given as a number or a list of numbers; specifies those boundaries where weak conditions will be enforced
- `pde`:
  - `formulation` - `default` or `lavrentovich`
  - `grad_desc` - If true, will use gradient descent to solve PDE
  - `gd_ls` - `backtrack`, `exact1`, `exact2`, or `none` - Chooses the line search type of the gradient descent
  - `solver` - `dynamic`, `newton`, or `builtin_nonlinear`
  - `tol_u` - lower tolerance for dynamic solver
  - `tol_l` - upper tolerance for dynamic solver
  - `tol` - tolerance for finding a solution
- `solver`:
  - `ksp_type` - Krylov subspace method
  - `pc_type` - Preconditioner type
- `time`:
  - `num` - Number of time steps
  - `save_every` - Number of time steps at which to save the simulation
  - `step` - The step size
- `vis`:
  - `normal` - `upward` or `outwrd` - which direction the normal vector is pointing

## Choosing a mesh

Will fill in later.