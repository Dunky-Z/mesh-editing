#include "mesh_editing.h"

MeshEdit::MeshEdit()
{
}

MeshEdit::~MeshEdit()
{
}


void MeshEdit::CaculateLaplacianCotMatrix(const SurfaceMesh& mesh)
{
	std::vector<T> tripletlist;
	tripletlist.reserve(20);
	const int p_num = mesh.n_vertices();
	L.resize(p_num, p_num);
	for (auto fit : mesh.faces())
	{
		Point p[3];
		double cot[3];
		int id[3];
		auto vf  = mesh.vertices(fit);
		for(int i = 0; i < 3; ++i, ++vf)
		{
			p[i] = mesh.position(*vf);
			id[i] = (*vf).idx();
		}
		for (int i = 0; i < 3; ++i)
		{
			int j = (i + 1) % 3, k = (j + 1) % 3;
			cot[i] = dot(p[j] - p[i], p[k] - p[i]) / norm(cross(p[j] - p[i], p[k] - p[i]));

			tripletlist.push_back({ id[j], id[k], -0.5 * cot[i] });
			tripletlist.push_back({ id[k], id[j], -0.5 * cot[i] });
		}
		for (int i = 0; i < 3; ++i)
		{
			tripletlist.push_back({ id[i], id[i], 0.5*(cot[(i + 1) % 3] + cot[(i + 2) % 3]) });
		}
	}
	L.setFromTriplets(tripletlist.begin(), tripletlist.end());
}

void MeshEdit::GetGradientField(const SurfaceMesh & mesh)
{

}
