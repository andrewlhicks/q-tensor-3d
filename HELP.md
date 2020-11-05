# Help

## Settings configuration

All settings files are contained in the folder 'settings'. Settings files are in the '.yml' format. Visit [PyYAML](https://pyyaml.org/wiki/PyYAMLDocumentation) for more information.

- `const` - The constants that will be used in the PDE itself
  - `L0` - Convex splitting constant
  - `L1` - Elastic energy constant L1
  - `L2` - Elastic energy constant L2
  - `L3` - Elastic energy constant L3
  - `A` - Bulk energy constant A
  - `B` - Bulk energy constant B
  - `C` - Bulk energy constant C
  - `ep` - The bulk energy is scaled by a factor of 1/ep^2. So as ep gets smaller, the end time needs to be increased

- `meshdata` - Specific data corresponding to the mesh being used
  - `file_path` - Specify the path to the mesh .msh file
  - `numnodes_init` - For manufactured solution, the number of nodes in each dimension of the unit cube mesh
  - `numnodes_max` - For manufactured solution, the maximum allowable number of nodes in each dimension of the unit cube mesh

- `options` - Take boolean values
  - `omit_init_printoff` - If 'true', omits the initial printoff of the settings
  - `visualize` - If 'true', creates a Paraview file to visualize the data
  - `manufactured` - If 'true', manufactures a solution, creates a unit cube mesh, and loops through different numbers of degrees of freedom

- `visdata` - Specifies any paraview settings that need to be configured
  - `file_path` - paraview/q-tensor-3d.pvd

- `solverdata` - Settings specific to the PDE solver being used
  - `ksp_type` - Krylov subspace method, use 'cg' (conjugate gradient) for symmetric positive definite matrices
  - `pc_type` - Preconditioner type

- `timedata` - Specifies the time step of the time-dependent PDE and the end time
  - `time_step`
  - `end_time`

## Choosing a mesh

If using a manufactured solution, a mesh covering the unit cube, i.e. the cube with diagonal from (0,0,0) to (1,1,1) in the Cartesian coordinate space, will be used automatically. The boundary edges in this mesh are numbered as follows:

1. The plane x = -1
2. The plane x = 1
3. The plane y = -1
4. The plane y = 1

This mesh will have the same number of nodes in each dimension. The number of nodes can be changed in settings by modifying `meshdata.numnodes_init` and `meshdata.numnodes_max`. For example, if the number of nodes is set to 10, then the resulting mesh will have 10 x 10 x 10 = 1000 nodes.

If a manufactured solution isn't being used, then you need to use a custom mesh file, which should be in the .msh format. Specify the path to this mesh in `settings.meshdata.file_path`.