
#include "mesh_editing.h"

int main()
{
	MeshEdit meshedit;
	SurfaceMesh mesh;
	mesh.read(meshedit.ori_mesh_path);
	cout << mesh.n_vertices() << endl;



	getchar();
	return 0;
}
