#pragma once
#include <iostream>

#include "pmp/SurfaceMesh.h"
#include "pmp/algorithms/DifferentialGeometry.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Cholesky>

using namespace Eigen;
using namespace pmp;
using namespace std;
using Point = pmp::Point;

typedef Eigen::Triplet<double> T;


class MeshEdit
{
public:
	MeshEdit();
	~MeshEdit();
	void CaculateLaplacianCotMatrix(const SurfaceMesh& mesh, Eigen::SparseMatrix<double> & L);
	pmp::vec3 MeshEdit::CaculateGradientField(const SurfaceMesh & mesh, const Face & face);
	pmp::vec3 CaculateGradientBu(const SurfaceMesh & mesh, const Face &face, const Vertex & vertex);
	Eigen::MatrixXf CaculateDivergence(SurfaceMesh & mesh, const Face & face);
	void SolvePoissonFunction(SurfaceMesh & mesh, const Eigen::SparseMatrix<double> & A, const Eigen::MatrixXf&  b);
	void Apply();
public:
	const string res_mesh_path = "D:/ITabc/ITabc/mesh-editing/res/res_mesh01.obj";
	const string ori_mesh_path = "D:/ITabc/ITabc/mesh-editing/build/model/ori.obj";
private:

};

