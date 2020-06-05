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
*/
void MeshEdit::CaculateLaplacianCotMatrix(const SurfaceMesh& mesh, Eigen::SparseMatrix<float> & L)
{
	std::vector<Tri> tripletlist;
	tripletlist.reserve(20);
	const int p_num = mesh.n_vertices();
	L.resize(p_num, p_num);
	for (auto fit : mesh.faces())
	{
		vec3 p[3];
		float cot[3];
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

			tripletlist.push_back(Tri(id[j], id[k], -0.5 * cot[i]));
			tripletlist.push_back(Tri(id[k], id[j], -0.5 * cot[i]));
		}
		for (int i = 0; i < 3; ++i)
		{
			tripletlist.push_back(Tri(id[i], id[i], 0.5*(cot[(i + 1) % 3] + cot[(i + 2) % 3])));
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
*/
pmp::vec3 MeshEdit::CaculateGradientBu(const SurfaceMesh & mesh, const Face &face, const Vertex & vertex)
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

	//\nabla{f(u)} = f_i\nabla{B_i(u)} + f_j\nabla{B_j(u)} + f_k\nabla{B_k(u)}
	fu[0] = points[0][0] * CaculateGradientBu(mesh, face, vertex[0])[0] + points[1][0] * CaculateGradientBu(mesh, face, vertex[1])[0] + points[2][0] * CaculateGradientBu(mesh, face, vertex[2])[0];
	fu[1] = points[0][1] * CaculateGradientBu(mesh, face, vertex[0])[1] + points[1][1] * CaculateGradientBu(mesh, face, vertex[1])[1] + points[2][1] * CaculateGradientBu(mesh, face, vertex[2])[1];
	fu[2] = points[0][2] * CaculateGradientBu(mesh, face, vertex[0])[2] + points[1][2] * CaculateGradientBu(mesh, face, vertex[1])[2] + points[2][2] * CaculateGradientBu(mesh, face, vertex[2])[2];

	return fu;
}


/*!
*@brief  求每个顶点的散度
*@param[out]
*@param[in]  SurfaceMesh & mesh	
*@return     Eigen::MatrixX3d	顶点的三个方向的散度
*/Eigen::MatrixXf MeshEdit::CaculateDivergence(const SurfaceMesh & mesh)
{
	int num_constrain = constraint_idx.size(), num_move = move_idx.size();
	size_t num_vertex = mesh.n_vertices();
	Eigen::MatrixXf div_w(num_vertex, 3);

	for (auto v : mesh.vertices())
	{
		vec3 div = vec3(0.0, 0.0, 0.0);
		for (auto fv : SurfaceMesh::FaceAroundVertexCirculator(&mesh, v))
		{
			float area = triangle_area(mesh, fv);
			vec3 bu = CaculateGradientBu(mesh, fv, v);
			vec3 grad = CaculateGradientField(mesh, fv);
			for (int i = 0; i < 3; ++i)
			{
				div[i] += area * bu[i] * grad[i];
			}
		}
		div_w.coeffRef(v.idx(), 0) = div[0];
		div_w.coeffRef(v.idx(), 1) = div[1];
		div_w.coeffRef(v.idx(), 2) = div[2];
	}

	div_w.conservativeResize(num_vertex + num_constrain + num_move, 3);
	std::cout << div_w.rows() << std::endl;

	// 用形变前坐标对固定锚点坐标进行赋值
	for (auto i = 0; i < num_constrain; i++)
	{
		div_w.coeffRef(i + num_vertex, 0) = constrain_pos[i][0];
		div_w.coeffRef(i + num_vertex, 1) = constrain_pos[i][1];
		div_w.coeffRef(i + num_vertex, 2) = constrain_pos[i][2];
	}
	int expand = 0.1;
	// 用形变后坐标对移动锚点坐标进行赋值
	for (auto i = 0; i < num_move; i++)
	{
		div_w.coeffRef(i + num_vertex + num_constrain, 0) = move_pos[i][0];
		div_w.coeffRef(i + num_vertex + num_constrain, 1) = move_pos[i][1];
		div_w.coeffRef(i + num_vertex + num_constrain, 2) = move_pos[i][2] + expand;
	}

	return div_w;
}

/*!
*@brief
*@param[out]
*@param[in]  SurfaceMesh & mesh
*@param[in]  const Eigen::SparseMatrix<float> & A
*@param[in]  const Eigen::MatrixXf & b
*@return     void
*/
void MeshEdit::SolvePoissonFunction(SurfaceMesh & mesh, const Eigen::SparseMatrix<float> & L, Eigen::MatrixXf&  b)
{
	Eigen::SparseMatrix<float> A = L;
	int num_vertex = mesh.n_vertices();
	int num_constrain = constraint_idx.size(), num_move = move_idx.size();
	//将A矩阵扩展而保持原有数据不变
	A.conservativeResize(num_vertex + num_constrain + num_move, num_vertex);
	cout << A.rows() << endl;
	cout << A.cols() << endl;

	for (auto i = 0; i < num_constrain; i++)
	{
		for (auto j = 0; j < num_vertex; j++)
		{
			if (j == constraint_idx[i])
				A.coeffRef(num_vertex + i, j) = 1;
			else
				A.coeffRef(num_vertex + i, j) = 0;
		}
	}

	// 移动锚点
	for (auto i = 0; i < num_move; i++)
	{
		for (auto j = 0; j < num_vertex; j++)
		{
			if (j == move_idx[i])
				A.coeffRef(num_vertex + num_constrain + i, j) = 1;
			else
				A.coeffRef(num_vertex + num_constrain + i, j) = 0;
			//cout << j << endl;
		}
	}

	A.makeCompressed();
	Eigen::SparseMatrix<float> AT = A.transpose();
	Eigen::SparseMatrix<float> ATA = AT * A;
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<float> > solver;
	solver.compute(ATA);
	Eigen::MatrixXf new_vertice(mesh.n_vertices(), 3);
	new_vertice.col(0) = solver.solve(AT*b.col(0));
	new_vertice.col(1) = solver.solve(AT*b.col(1));
	new_vertice.col(2) = solver.solve(AT*b.col(2));

	for (auto v : mesh.vertices())
	{
		Point& pos = mesh.position(v);
		pos[0] = new_vertice.coeff(v.idx(), 0);
		pos[1] = new_vertice.coeff(v.idx(), 1);
		pos[2] = new_vertice.coeff(v.idx(), 2);
	}
	mesh.write(res_mesh_path);
}


/*!
*@brief  设置约束的固定顶点和移动顶点
*@param[out]
*@param[in]  SurfaceMesh & mesh
*@return     void
*/void MeshEdit::SetConstraintAndMoveVertex(const SurfaceMesh & mesh)
{
	for (auto v : mesh.vertices())
	{
		//--------------cubic---------------///
		////设置z坐标小于-0.98的顶点为约束顶点。大概就是脚部
		//if (mesh.position(v)[2] == 1.0)
		//{
		//	constraint_idx.push_back(v.idx());
		//	constrain_pos.push_back(mesh.position(v));
		//}
		////设置z坐标大于0.63的顶点为可移动顶点，大概是头部位置
		//if (mesh.position(v)[2] == -1.0)
		//{
		//	move_idx.push_back(v.idx());
		//	move_pos.push_back(mesh.position(v));
		//}

		//--------------Cylinder---------------///
		////圆柱体z坐标小于11固定
		//if (mesh.position(v)[2] < 11.0)
		//{
		//	constraint_idx.push_back(v.idx());
		//	constrain_pos.push_back(mesh.position(v));
		//}
		////设置z坐标大于18为移动
		//if (mesh.position(v)[2] > 18.0)
		//{
		//	move_idx.push_back(v.idx());
		//	move_pos.push_back(mesh.position(v));
		//}

		//--------------ori---------------///
		//设置z坐标小于-0.98的顶点为约束顶点。大概就是脚部
		if (mesh.position(v)[2] < -0.97)
		{
			constraint_idx.push_back(v.idx());
			constrain_pos.push_back(mesh.position(v));
		}
		//设置z坐标大于0.63的顶点为可移动顶点，大概是头部位置
		if (mesh.position(v)[2] > 0.6)
		{
			move_idx.push_back(v.idx());
			move_pos.push_back(mesh.position(v));
		}
	}
}


/*!
*@brief  主函数
*@param[out] 
*@return     void  
*/void MeshEdit::Apply()
{
	MeshEdit meshedit;
	SurfaceMesh mesh;
	mesh.read(ori_mesh_path);
	cout << mesh.n_vertices() << endl;
	int num_vertex = mesh.n_vertices();
	Eigen::SparseMatrix<float> A(num_vertex, num_vertex);

	meshedit.SetConstraintAndMoveVertex(mesh);
	Eigen::MatrixXf b = meshedit.CaculateDivergence(mesh);
	meshedit.CaculateLaplacianCotMatrix(mesh, A);
	meshedit.SolvePoissonFunction(mesh, A, b);
}

