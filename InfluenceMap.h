#pragma once

#include <vector>
#include <Eigen/Geometry>
#include <string>

void write_influence_color_map_OBJ(const std::string& file_name, const Eigen::MatrixXd& V,
	const Eigen::MatrixXi& T, const Eigen::MatrixXd& W, const std::vector<int>& control_vertices_idx, int cage_vertices_offset, bool transposeW);