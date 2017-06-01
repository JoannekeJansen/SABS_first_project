from dolfin import *

#mesh = RectangleMesh(Point(0.0, 0.0), Point(12.0, 0.01), 200, 1)
mesh = UnitSquareMesh(2, 2)

gdim = mesh.geometry().dim()
Q = FunctionSpace(mesh, "CG", 1)
dofmap = Q.dofmap()
dofs = dofmap.dofs()

f = Expression('x[0]*x[1]',degree=2)

f_proj = project(f, Q)

F = f_proj.vector().array()

# Get coordinates as len(dofs) x gdim array
dofs_x = Q.tabulate_dof_coordinates().reshape((-1, gdim))

for dof, dof_x in zip(dofs, dofs_x):
    print dof, ':', dof_x
    print dof_x[0]
    print F[dof]

#print 'dof index | dof coordinate |  dof value'
for i, (x, v) in enumerate(zip(dofs_x, F)):
    print i, x, v 
