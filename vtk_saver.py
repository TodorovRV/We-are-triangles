import vtk
from matplotlib import pyplot as plt


def make_snapshot(points, triangles, j):
    unstructuredGrid = vtk.vtkUnstructuredGrid()
    vtk_points = vtk.vtkPoints()
    for i in range(int(len(points) / 3)):
        x = points[i * 3]
        y = points[i * 3 + 1]
        z = points[i * 3 + 2]
        vtk_points.InsertNextPoint(x, y, z)
        unstructuredGrid.SetPoints(vtk_points)
    for i in range(int(len(triangles) / 3)):
        tr = vtk.vtkTriangle()
        tr.GetPointIds().SetId(0, triangles[i * 3])
        tr.GetPointIds().SetId(1, triangles[i * 3 + 1])
        tr.GetPointIds().SetId(2, triangles[i * 3 + 2])
        unstructuredGrid.InsertNextCell(tr.GetCellType(), tr.GetPointIds())
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetInputDataObject(unstructuredGrid)
    writer.SetFileName("Snapshots/surface" + str(j) + ".vtu")
    writer.Write()


points = [0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0]
triangles = [0, 1, 2, 2, 1, 3]
calc_points = []
with open('output.txt') as f:
    calc_triangles = list(map(int, f.readline().split()))
    iterations = int(f.readline())
    for i in range(iterations):
        calc_points.append(list(map(float, f.readline().split())))
# x = []
# y = []
# for i in range(int(len(calc_points) / 3)):
#     x.append(calc_points[i * 3])
#     y.append(calc_points[i * 3 + 1])
# plt.scatter(x, y)
# plt.show()
print(len(calc_points))
for i in range(len(calc_points)):
    make_snapshot(calc_points[i], calc_triangles, i)
#make_snapshot(points, triangles)






