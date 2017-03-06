std::vector<Eigen::Matrix3d> calculateS(const Eigen::MatrixXd &vertices, const Eigen::MatrixXi &faces)
{
	std::vector<Eigen::Matrix3d> S;

	for (size_t i = 0, rows = faces.rows(); i < rows; ++i)
	{
		Eigen::Vector3d v1, v2, v3, v4;
		
		v1(0) = vertices(faces(i, 0), 0);
		v1(1) = vertices(faces(i, 0), 1);
		v1(2) = vertices(faces(i, 0), 2);
		
		v2(0) = vertices(faces(i, 1), 0);
		v2(1) = vertices(faces(i, 1), 1);
		v2(2) = vertices(faces(i, 1), 2);
		
		v3(0) = vertices(faces(i, 2), 0);
		v3(1) = vertices(faces(i, 2), 1);
		v3(2) = vertices(faces(i, 2), 2);
		
		v4 = (v2 - v1).cross(v3 - v1);
		v4 = v1 + v4 / sqrt(v4.norm());
		
		Eigen::Matrix3d s;
		s.col(0) = v2 - v1;
		s.col(1) = v3 - v1;
		s.col(2) = v4 - v1;
		
		S.push_back(s);
	}
	
	return S;
}


std::vector<Eigen::Matrix3d> calculateTransformation(const std::vector<Eigen::Matrix3d> &S0, const std::vector<Eigen::Matrix3d> &Si)
{
	std::vector<Eigen::Matrix3d> S;
	
	for (size_t i = 0, tetrahedra = S0.size(); i < tetrahedra; ++i)
	{
		S.push_back(Si[i] * S0[i].inverse());
	}

	return S;
}


Eigen::MatrixXd buildG(const size_t number_of_triangles)
{
	size_t rows = 6 * number_of_triangles;
	Eigen::MatrixXd G = Eigen::MatrixXd::Zero(rows, 9*number_of_triangles);

	size_t i = 0, j = 0;
	while (i < rows)
	{
		G(i, j) = -1;
		G(i, j+3) = 1;
		++i;
		G(i, j+1) = -1;
		G(i, j+4) = 1;
		++i;
		G(i, j+2) = -1;
		G(i, j+5) = 1;
		++i;
		G(i, j) = -1;
		G(i, j+6) = 1;
		++i;
		G(i, j+1) = -1;
		G(i, j+7) = 1;
		++i;
		G(i, j+2) = -1;
		G(i, j+8) = 1;
		++i;

		j += 9;
	}

	return G;
}


Eigen::MatrixXd buildH(const std::vector<Eigen::Matrix3d> S)
{
	size_t dimension = 6 * S.size();
	Eigen::MatrixXd H = Eigen::MatrixXd::Zero(dimension, dimension);

	size_t i = 0, j = 0, k = 0;
	while (i < dimension)
	{
		for (size_t l = 0; l < 3; ++l)
		{
			H(i, j) = S[k](l, 0);
			H(i, j+1) = S[k](l, 1);
			H(i, j+2) = S[k](l, 2);
			++i;
		}

		j += 3;
		for (size_t l = 0; l < 3; ++l)
		{
			H(i, j) = S[k](l, 0);
			H(i, j+1) = S[k](l, 1);
			H(i, j+2) = S[k](l, 2);
			++i;
		}

		j += 3;
		++k;
	}

	return H;
}


Eigen::VectorXd buildB(const Eigen::MatrixXd &vertices, const Eigen::MatrixXi &faces)
{
	size_t number_of_triangles = faces.rows();
	Eigen::VectorXd b(9*number_of_triangles);

	size_t i = 0;
	for (size_t j = 0; j < number_of_triangles; ++j)
	{
		for (size_t k = 0; k < 3; ++k)
		{
			b(i++) = vertices(faces(j, k), 0);
			b(i++) = vertices(faces(j, k), 1);
			b(i++) = vertices(faces(j, k), 2);
		}
	}

	return b;
}


Eigen::MatrixXd pinv(const Eigen::MatrixXd &m)
{
	size_t rows = m.rows(), cols = m.cols();

	arma::mat M(rows, cols);
	for (size_t i = 0; i < rows; ++i)
		for (size_t j = 0; j < cols; ++j)
			M(i, j) = m(i, j);

	arma::mat pseudoinverse = arma::pinv(M);
	rows = pseudoinverse.n_rows;
	cols = pseudoinverse.n_cols;

	Eigen::MatrixXd pinv(rows, cols);
	for (size_t i = 0; i < rows; ++i)
		for (size_t j = 0; j < cols; ++j)
			pinv(i, j) = pseudoinverse(i, j);

	return pinv;
}


Eigen::MatrixXd bToVertices(const Eigen::VectorXd &b, const Eigen::MatrixXi &faces,
	const size_t number_of_vertices)
{
	Eigen::MatrixXd vertices(number_of_vertices, 3);
	
	size_t j = 0;
	for (size_t i = 0, rows = faces.rows(); i < rows; ++i)
	{
		vertices(faces(i, 0), 0) = b(j++);
		vertices(faces(i, 0), 1) = b(j++);
		vertices(faces(i, 0), 2) = b(j++);

		vertices(faces(i, 1), 0) = b(j++);
		vertices(faces(i, 1), 1) = b(j++);
		vertices(faces(i, 1), 2) = b(j++);

		vertices(faces(i, 2), 0) = b(j++);
		vertices(faces(i, 2), 1) = b(j++);
		vertices(faces(i, 2), 2) = b(j++);
	}

	return vertices;
}
