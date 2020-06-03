#include "mesh_editing.h"

MeshEdit::MeshEdit()
{
}

MeshEdit::~MeshEdit()
{
}


/*!
*@brief  计算拉普拉斯权值矩阵
*@param[out] 拉普拉斯稀疏矩阵
*@param[in]  const SurfaceMesh & mesh  待求拉普拉斯矩阵的原始网格
*@param[in]  Eigen::SparseMatrix<float> & L  拉普拉斯矩阵，是个大型稀疏矩阵
*@return     void
*/void MeshEdit::CaculateLaplacianCotMatrix(const SurfaceMesh& mesh, Eigen::SparseMatrix<double> & L)
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
		auto vf = mesh.vertices(fit);
		for (int i = 0; i < 3; ++i, ++vf)
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


/*!
*@brief  计算基函数B(u)的梯度。输入面片和顶点信息，在该面片上计算出输入顶点所对应的基函数梯度信息。
*@param[out]
*@param[in]  const SurfaceMesh & mesh  原始网格信息
*@param[in]  const Face & face  某一个待求解的面片
*@param[in]  const Vertex & vertex  面片face上的一个顶点
*@return     pmp::vec3  B(u)的梯度
*/pmp::vec3 MeshEdit::CaculateGradientBu(const SurfaceMesh & mesh, const Face &face, const Vertex & vertex)
{
	std::vector<Point> points;
	int idx = 0;
	int i = 0;
	for (auto vf : mesh.vertices(face))
	{
		if (vf == vertex)
		{
			idx = i;
		}
		points.push_back(mesh.position(vf));
		++i;
	}

	vec3 face_normal = cross(points[1] - points[0], points[2] - points[0]);
	vec3 face_normal_normalized = normalize(face_normal);
	float face_area = 0.5 * norm(face_normal);
	//和顶点vertex相对的那条边
	vec3 edge_v_opposite;
	edge_v_opposite = points[(idx + 1) % 3] - points[(idx + 2) % 3];
	//通过面法向量和对边向量叉乘求对边向量旋转九十度的向量
	vec3 gradient_direction = cross(face_normal_normalized, edge_v_opposite);
	vec3 bu = gradient_direction / (2 * face_area);
	return bu;
}


/*!
*@brief  计算面片的梯度场
*@param[out]
*@param[in]  const SurfaceMesh & mesh  输入源网格
*@param[in]  const Face & face
*@return     pmp::vec3
*/
pmp::vec3 MeshEdit::CaculateGradientField(const SurfaceMesh & mesh, const Face & face)
{
	std::vector<Point> points;
	std::vector<Vertex> vertex;
	for (auto vf : mesh.vertices(face))
	{
		points.push_back(mesh.position(vf));
		vertex.push_back(vf);
	}
	vec3 fu;

	//fb_j += (points[1][i] - points[0][i])*CaculateGradientBu(mesh, face, vertex[1]);
	//fb_k += (points[2][i] - points[0][i])*CaculateGradientBu(mesh, face, vertex[2]);
	fu[0] = points[0][0] * CaculateGradientBu(mesh, face, vertex[0])[0] + points[1][0] * CaculateGradientBu(mesh, face, vertex[0])[0]+ points[2][0] * CaculateGradientBu(mesh, face, vertex[0])[0];
	fu[1] = points[0][1] * CaculateGradientBu(mesh, face, vertex[1])[0] + points[1][1] * CaculateGradientBu(mesh, face, vertex[1])[0]+ points[2][1] * CaculateGradientBu(mesh, face, vertex[1])[1];
	fu[2] = points[0][2] * CaculateGradientBu(mesh, face, vertex[2])[0] + points[1][2] * CaculateGradientBu(mesh, face, vertex[2])[0]+ points[2][2] * CaculateGradientBu(mesh, face, vertex[2])[2];


	//(fj-fi)*bu_j+(fk-fi)*bu_k
	return fu;
}


Eigen::MatrixXf MeshEdit::CaculateDivergence(SurfaceMesh & mesh, const Face & face)
{
	Eigen::VectorXf div_w;
	size_t n = mesh.n_vertices();
	div_w.setZero(n);
	int index = 0;
	int id = 0;
	for (auto v : mesh.vertices())
	{
		if (id == 5359)
		{
			Point &pos = mesh.position(v);
			pos = Point{ 0.326361, -0.278979, -0.140672 };
		}
		id++;

	}

	for (auto v : mesh.vertices())
	{
		Scalar div = 0;
		for (auto fv : SurfaceMesh::FaceAroundVertexCirculator(&mesh, v))
		{
			Scalar area = triangle_area(mesh, fv);
			vec3 bu = CaculateGradientBu(mesh, fv, v);
			for (int i = 0; i < 3; ++i)
			{
				div += area * bu[i] * CaculateGradientField(mesh, fv)[i];
			}
		}
		div_w[index++] = div;
	}
	return div_w;
}

void MeshEdit::SolvePoissonFunction(SurfaceMesh & mesh, const Eigen::SparseMatrix<double> & A, const Eigen::MatrixXf&  b)
{
	Eigen::SparseMatrix<double> AT = A.transpose();
	Eigen::SparseMatrix<double> ATA = AT * A;
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double> > MatricsCholesky(ATA);
	Eigen::VectorXf new_vertice = MatricsCholesky.solve(AT * b);

	mesh.write(res_mesh_path);
}


void MeshEdit::Apply()
{
	MeshEdit meshedit;
	SurfaceMesh mesh;
	mesh.read(ori_mesh_path);
	cout << mesh.n_vertices() << endl;

	Eigen::SparseMatrix<float> A;

	//CaculateLaplacianCotMatrix(mesh, A);

	//cout << A.topLeftCorner(50, 50) << endl;
	for (auto v : mesh.vertices())
	{
		for (auto f : mesh.faces(v))
		{
			vec3 bu = meshedit.CaculateGradientBu(mesh, f, v);
			cout << "bu :         " << bu << endl;
			vec3 fu = meshedit.CaculateGradientField(mesh, f);
			cout << "fu:    " << fu << endl;
			Eigen::MatrixXf dw = CaculateDivergence(mesh, f);
			cout << dw.topRows(20)<< endl;
		}
	}
}