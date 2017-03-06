#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <armadillo>
#include "t-star.hpp"



int main()
{
	Eigen::MatrixXd vertices_source_neutral, vertices_source_smile, vertices_target_neutral;
	Eigen::MatrixXi faces_source_neutral, faces_source_smile, faces_target_neutral;
	igl::readOBJ("./test/neutral.obj", vertices_source_neutral, faces_source_neutral);
	igl::readOBJ("./test/smile.obj", vertices_source_smile, faces_source_smile);
	igl::readOBJ("./test/ugly.obj", vertices_target_neutral, faces_target_neutral);

	std::vector<Eigen::Matrix3d> S_source_neutral, S_source_smile, S_source;
	S_source_neutral = calculateS(vertices_source_neutral, faces_source_neutral);
	S_source_smile = calculateS(vertices_source_smile, faces_source_smile);
	S_source = calculateTransformation(S_source_neutral, S_source_smile);
	
	Eigen::MatrixXd H = buildH(S_source);
	Eigen::MatrixXd G = buildG(faces_source_neutral.rows());
	Eigen::MatrixXd G_pinv = pinv(G);
	Eigen::MatrixXd T_star = G_pinv * H * G;

	Eigen::VectorXd b0_source = buildB(vertices_source_neutral, faces_source_neutral);
	Eigen::VectorXd bi_source = buildB(vertices_source_smile, faces_source_smile);
	Eigen::VectorXd b0_target = buildB(vertices_target_neutral, faces_target_neutral);

	Eigen::VectorXd translation = b0_target - b0_source + bi_source;
	translation = translation - G_pinv * G * translation;
	
	Eigen::VectorXd bi_target = T_star * b0_target + translation;
	
	Eigen::MatrixXd vertices_target_smile = bToVertices(
													bi_target, faces_target_neutral, vertices_target_neutral.rows());
	igl::writeOBJ("./test/output_ugly.obj", vertices_target_smile, faces_target_neutral);

	return 0;
}
