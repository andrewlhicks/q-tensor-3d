# Q-tensor 3D

"Q-tensor 3D" is an implementation of the Landau-de Gennes Q-tensor model for liquid crystals.

## Saves

To use qtensor3d, a save folder needs to be created. The easiest way to do this is to use `qsave` in the following manner:
```
qsave -b <savepath>
```
where `<savepath>` is the save folder. This will create a save folder with the default settings, cosntants, and uflexpr files.

## Settings configuration

Settings files are in the '.yml' format. Visit [PyYAML](https://pyyaml.org/wiki/PyYAMLDocumentation) for more information.

- `mesh`:
  - `source` - `builtin`, `local`, `global`, or `legacy`
  - `name` - See "Choosing a mesh" for more info
  - `refs` - How many refinements of the mesh to complete
- `options`:
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

To choose a mesh, we go into our `settings.yml` file and edit the settings listed under `mesh`.
We begin with the `source` setting.
We may choose `builtin`, `local`, `global`, or `legacy`.

### Builtin

If using a builtin mesh, the only mesh implemented at the moment is the Box Mesh.
Let's say we want to create a Box Mesh of length, width and height `x`, `y`, and `z` respectively.
Moreover, let's let `xn`, `yn`, and `zn` be the number of nodes along the length, width, and height repectively.
To do this, we set `source` to `builtin`, and set `name` to `BoxMesh xn yn zn x y z`.
For example, for a slab with a 10 x 10 base and a height of 0.2, with a mesh size of 0.1, we would implement it as follows:
```
mesh:
  source: builtin
  name: BoxMesh 10 10 2 1 1 0.2
```

### Local

To use a local mesh, simply save a `.msh` file inside of the save folder.
For example, if we want to use a mesh named `mesh.msh` located in our save folder, we implement the following settings:
```
mesh:
  source: local
  name: mesh.msh
```

### Global

To use a mesh saved at an absolute file path, specify the absolute path.
For example, if we have a mesh located at `/scratch/username/meshes/mesh.msh`, we would implement is as follows:
```
mesh:
  source: global
  name: /scratch/username/meshes/mesh.msh
```

### Legacy

Mainly used for the old mesh format. For instance, we would specify `/scratch/username/meshes/meshX.msh`, where `X` is the refinement number.
For example, the following implementation will cycle through 3 refinement levels of the mesh named as above:
```
mesh:
  source: legacy
  name: /scratch/username/meshes/mesh
  refs: 0-2
```