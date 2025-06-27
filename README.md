# Q-tensor 3D

"Q-tensor 3D" is an implementation of the Landau-de Gennes Q-tensor model for liquid crystals.

## Saves

To use qtensor3d, a save folder needs to be created. The easiest way to do this is to use `qsave` in the following manner:
```
qsave -b <savepath>
```
where `<savepath>` is the save folder. This will create a save folder with the default settings, cosntants, and uflexpr files.

## Understanding the settings file

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

## Understanding the User Expression file

The file `userexpr.yml` allows the user to input custom expressions for the initial condition _q-vector_ (`initcond`), the weak boundary _director_ (`w_bdy_nu`), the strong boundary _q-vector_ (`sbdy`), the manufactured _q-vector_ (`manu_q`), the forcing right hand side _q-vector_ on the bulk (`forcing_f`), and the forcing right hand side _q-vector_ on the boundary (`forcing_g`).

### Creating user expressions using UFL objects

Using UFL objects is now the preferred method to specify user expressions.
The two major classes of UFL objects that the user can specify are _q-vectors_ (using the `!qvector` flag) and _directors_ (using the `!director` flag).
Vector objects (3-dimensional) are only used for the weak boundary director (`w_bdy_nu`), while q-vector objects (5-dimensional) are used for everything else.

Besides the hundreds of already available objects in the standard UFL library, there are several additional functions available to the user, most of which can be nested inside of each other:

#### `as_vector`

Outputs a vector of any dimension; thus can be used to specify a q-vector or a director object.

__Example:__
```
initcond: !qvector as_vector([x0**2, x0+x1, x1+x2, 1, 0])
```
Another example:
```
w_bdy_nu: !director as_vector([0,0,1])
```

#### `from_director`

Outputs a q-vector from a given director.

__Example:__
```
initcond: !qvector from_director([cos(5*x2),sin(5*x2),0])
```

#### `from_spherical_director`

Outputs a q-vector from a given director with _spherical_ coordinates.
Equivalent to passing the director through `spherical_to_cartesian` and then to `from_director`.

Example:
```
initcond: !qvector from_spherical_director([0,1,0])
```

#### `spherical_to_cartesian`

This converts a 3-dimensional spherical vector to the corresponding 3-dimensional Cartesian vector.
It is written as
```
spherical_to_cartesian([a,b,c])
```
where `a`, `b`, and `c` are arbitrary scalar expressions.
From here a vector `v` is created:
```
v  =  a * e_r  +  b * e_theta  +  c * e_phi
```
where
```
e_r = as_vector([x0, x1, x2])/sqrt(x0**2+x1**2+x2**2)
e_theta = as_vector([x0*x2, x1*x2, -(x0**2+x1**2)])/sqrt((x0**2+x1**2+x2**2)*(x0**2+x1**2))
e_phi = as_vector([-x1, x0, 0])/sqrt(x0**2+x1**2)
```

### `smooth_transition`

Invoking the function `smooth_transition(x, I=[a,b])` will give a smooth transition (scalar) function in the variable `x` on the interval `I = [a,b]`, with the value of `0` at `x=a` and `1` at `x=b`.

__Example:__ Suppose one wanted to define a function that is equal to `0` at all values of `x2` less than or equal to `-2`, equal to `1` at all values of `x2` greater than or equal to `3`, and with a smooth transition between the two values on the interval `(-2,3)`.
In this case one would write
```
smooth_transition(x2, I=[-2,3])
```

### Creating user expressions using Sympy objects (deprecated)

Will add documentation later, though please note that this is now deprecated and the preferred method of creating user expressions is through UFL objects.