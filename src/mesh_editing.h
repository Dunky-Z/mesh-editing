#pragma once
#include <iostream>

#include "pmp/SurfaceMesh.h"
#include "pmp/algorithms/DifferentialGeometry.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Cholesky>

using namespace pmp;
using namespace std;
using namespace Eigen;


using Point = pmp::Point;
typedef Eigen::Triplet<float> Tri;


class MeshEdit
{
public:
	MeshEdit();
	~MeshEdit();
	void CaculateLaplacianCotMatrix(const SurfaceMesh& mesh, Eigen::SparseMatrix<float> & L);
	pmp::vec3 MeshEdit::CaculateGradientField(const SurfaceMesh & mesh, const Face & face);
	pmp::vec3 CaculateGradientBu(const SurfaceMesh & mesh, const Face &face, const Vertex & vertex);
	Eigen::MatrixXf CaculateDivergence(const SurfaceMesh & mesh);
	void SolvePoissonFunction(SurfaceMesh & mesh, const Eigen::SparseMatrix<float> & A, Eigen::MatrixXf&  b);
	void Apply();
	void SetConstraintAndMoveVertex(const SurfaceMesh & mesh);
public:
	const string res_mesh_path = "D:/ITabc/ITabc/mesh-editing/res/res_mesh03.obj";
	//const string ori_mesh_path = "D:/ITabc/ITabc/mesh-editing/build/model/ori-remesh.obj";
	//const string ori_mesh_path = "D:/ITabc/ITabc/mesh-editing/build/model/Cylinder.obj";
	const string ori_mesh_path = "D:/ITabc/ITabc/mesh-editing/build/model/cube.obj";
	std::vector<pmp::vec3> constrain_pos, move_pos;
	std::vector<int> constraint_idx, move_idx;
private:

};

