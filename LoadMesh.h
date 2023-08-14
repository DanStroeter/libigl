#pragma once 

#include <string>

#include <Eigen/Geometry>

bool load_mesh(std::string const & file_name, Eigen::MatrixXd& V, Eigen::MatrixXi & T, double scaling_factor);

bool load_cage(std::string const& file_name, Eigen::MatrixXd& V, Eigen::VectorXi& P, Eigen::MatrixXi& CF,
	double scaling_factor, bool triangulate_quads, Eigen::MatrixXd* V_embedding = nullptr);
