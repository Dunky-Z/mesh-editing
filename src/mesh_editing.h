#pragma once
#include <iostream>

#include "pmp/SurfaceMesh.h"

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
	void CaculateLaplacianCotMatrix(const SurfaceMesh& mesh);
	void GetGradientField(const SurfaceMesh& mesh);
private:

};

