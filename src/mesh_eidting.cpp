#include "mesh_editing.h"

MeshEdit::MeshEdit()
{
}

MeshEdit::~MeshEdit()
{
}


/*!
*@brief  ����������˹Ȩֵ����
*@param[out] ������˹ϡ�����
*@param[in]  const SurfaceMesh & mesh  ����������˹�����ԭʼ����
*@param[in]  Eigen::SparseMatrix<float> & L  ������˹�����Ǹ�����ϡ�����
*@return     void
*/

void MeshEdit::CaculateLaplacianCotMatrix(const SurfaceMesh& mesh, Eigen::SparseMatrix<float> & L)
{
	//int vn = mesh.n_vertices();
	//int count0 = 0;
	//vector<int> begin_N(vn);
	//for (auto v : mesh.vertices()) {
	//	begin_N[v.idx()] = count0;
	//	count0 += mesh.valence(v) + 1;
	//}
	//typedef Eigen::Triplet<Scalar> T;
	//vector<T> tripletList(count0);
	////n*n
	//for (auto v : mesh.vertices()) {
	//	int i = 0;
	//	float sum_weight = 0;
	//	for (auto nei_v : SurfaceMesh::VertexAroundVertexCirculator(&mesh, v)) {
	//		Edge e = mesh.find_edge(v, nei_v);
	//		Scalar weight = 0.5 * cotan_weight(mesh, e);
	//		tripletList[begin_N[v.idx()] + i + 1] = T(v.idx(), nei_v.idx(), -weight);
	//		i++;
	//		sum_weight += weight;
	//	}
	//	tripletList[begin_N[v.idx()]] = T(v.idx(), v.idx(), sum_weight);
	//}
	//L.setFromTriplets(tripletList.begin(), tripletList.end());

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
*@brief  ���������B(u)���ݶȡ�������Ƭ�Ͷ�����Ϣ���ڸ���Ƭ�ϼ�������붥������Ӧ�Ļ������ݶ���Ϣ��
*@param[out]
*@param[in]  const SurfaceMesh & mesh  ԭʼ������Ϣ
*@param[in]  const Face & face  ĳһ����������Ƭ
*@param[in]  const Vertex & vertex  ��Ƭface�ϵ�һ������
*@return     pmp::vec3  B(u)���ݶ�
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
	//�Ͷ���vertex��Ե�������
	vec3 edge_v_opposite;
	edge_v_opposite = points[(idx + 1) % 3] - points[(idx + 2) % 3];
	//ͨ���淨�����ͶԱ����������Ա�������ת��ʮ�ȵ�����
	vec3 gradient_direction = cross(face_normal_normalized, edge_v_opposite);
	vec3 bu = gradient_direction / (2 * face_area);
	return bu;
}


/*!
*@brief  ������Ƭ���ݶȳ�
*@param[out]
*@param[in]  const SurfaceMesh & mesh  ����Դ����
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
*@brief  ��ÿ�������ɢ��
*@param[out]
*@param[in]  SurfaceMesh & mesh
*@return     Eigen::MatrixX3d
*/Eigen::MatrixXf MeshEdit::CaculateDivergence(SurfaceMesh & mesh)
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

	//------

	//vertex# 734	position[1.395040 6.366990 18.991699]
	//vertex# 420	position[1.578030 - 6.223970 18.994101]
	//vertex# 581	position[8.027470 - 0.384723 18.994600]
	//vertex# 1107	position[-4.404760 - 0.984303 18.991899]
	//vertex# 1158	position[2.181440 - 1.105770 19.000099]
	//vertex# 180	position[0.123524 - 0.299188 19.000000]
	//vertex# 1368	position[1.458380 - 0.537347 10.012600]

	////���������̶����һ���ƶ���
	//std::vector<std::vector<float>> fixed_pos{ {1.395040, 6.366990, 18.991699},{1.578030, -6.223970, 18.994101},{8.027470, -0.384723, 18.994600},{-4.404760, -0.984303 ,18.991899},{2.181440, -1.105770, 19.000099}, {0.123524, -0.299188 ,19.000000} };
	//std::vector<std::vector<float>> move_pos{ {1.458380, -0.537347, 5.0} };
	// ���α�ǰ����Թ̶�ê��������и�ֵ
	for (auto i = 0; i < num_constrain; i++)
	{
		div_w.coeffRef(i + num_vertex, 0) = constrain_pos[i][0];
		div_w.coeffRef(i + num_vertex, 1) = constrain_pos[i][1];
		div_w.coeffRef(i + num_vertex, 2) = constrain_pos[i][2];
	}
	int expand = 0.1;
	// ���α��������ƶ�ê��������и�ֵ
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
	//��A������չ������ԭ�����ݲ���
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

	// �ƶ�ê��
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


void MeshEdit::Apply()
{
	MeshEdit meshedit;
	SurfaceMesh mesh;
	mesh.read(ori_mesh_path);
	cout << mesh.n_vertices() << endl;
	int num_vertex = mesh.n_vertices();
	Eigen::SparseMatrix<float> A(num_vertex, num_vertex);

	//CaculateLaplacianCotMatrix(mesh, A);

	//cout << A.topLeftCorner(50, 50) << endl;
	//for (auto v : mesh.vertices())
	//{
	//	for (auto f : mesh.faces(v))
	//	{
	//		vec3 bu = meshedit.CaculateGradientBu(mesh, f, v);
	//		cout << "bu :         " << bu << endl;
	//		vec3 fu = meshedit.CaculateGradientField(mesh, f);
	//		cout << "fu:    " << fu << endl;
	//		Eigen::MatrixXf dw = CaculateDivergence(mesh, f);
	//		cout << dw.topRows(20)<< endl;
	//	}
	//}
	meshedit.SetConstraintAndMoveVertex(mesh);
	Eigen::MatrixXf b = meshedit.CaculateDivergence(mesh);
	meshedit.CaculateLaplacianCotMatrix(mesh, A);
	meshedit.SolvePoissonFunction(mesh, A, b);
}


/*!
*@brief  ����Լ���Ĺ̶�������ƶ�����
*@param[out]
*@param[in]  SurfaceMesh & mesh
*@return     void
*/void MeshEdit::SetConstraintAndMoveVertex(SurfaceMesh & mesh)
{
	for (auto v : mesh.vertices())
	{
		//--------------cubic---------------///
		////����z����С��-0.98�Ķ���ΪԼ�����㡣��ž��ǽŲ�
		//if (mesh.position(v)[2] == 1.0)
		//{
		//	constraint_idx.push_back(v.idx());
		//	constrain_pos.push_back(mesh.position(v));
		//}
		////����z�������0.63�Ķ���Ϊ���ƶ����㣬�����ͷ��λ��
		//if (mesh.position(v)[2] == -1.0)
		//{
		//	move_idx.push_back(v.idx());
		//	move_pos.push_back(mesh.position(v));
		//}

		//--------------Cylinder---------------///
		////Բ����z����С��11�̶�
		//if (mesh.position(v)[2] < 11.0)
		//{
		//	constraint_idx.push_back(v.idx());
		//	constrain_pos.push_back(mesh.position(v));
		//}
		////����z�������18Ϊ�ƶ�
		//if (mesh.position(v)[2] > 18.0)
		//{
		//	move_idx.push_back(v.idx());
		//	move_pos.push_back(mesh.position(v));
		//}

		//--------------ori---------------///
		//����z����С��-0.98�Ķ���ΪԼ�����㡣��ž��ǽŲ�
		if (mesh.position(v)[2] < -0.97)
		{
			constraint_idx.push_back(v.idx());
			constrain_pos.push_back(mesh.position(v));
		}
		//����z�������0.63�Ķ���Ϊ���ƶ����㣬�����ͷ��λ��
		if (mesh.position(v)[2] > 0.6)
		{
			move_idx.push_back(v.idx());
			move_pos.push_back(mesh.position(v));
		}
	}
}