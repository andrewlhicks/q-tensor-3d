from firedrake.mesh import MeshGeometry

from dataclasses import dataclass

__all__ = ('FunctionSpaceData',)

@dataclass
class FunctionSpaceData:
    mesh: MeshGeometry
    family: str
    degree: int

    def __iter__(self):
        return iter((self.mesh, self.family, self.degree)) # atrocious solution, but astuple() calls for a deep copy of self.mesh, and deep copy doesn't seem to work on instances of MeshGeometry