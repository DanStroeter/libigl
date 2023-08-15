#include "LoadMesh.h"

#include <vector>
#include <iostream>

#include <igl/readOBJ.h>
#include <igl/readMSH.h>
#include <igl/EPS.h>

template<typename vec, typename T>
void scale_verts(std::vector<vec>& verts, T fac)
{
	if (verts.empty())
	{
		return;
	}

	vec min_vert = verts[0], max_vert = verts[0];

	for (unsigned int i = 0; i < verts.size(); ++i)
	{
		auto const vert = verts[i];
		for (unsigned k = 0; k < 3u; ++k)
		{
			min_vert[k] = std::min(min_vert[k], vert[k]);
			max_vert[k] = std::max(max_vert[k], vert[k]);
		}
	}

	const vec mid = { (min_vert[0] + max_vert[0]) * .5f, (min_vert[1] + max_vert[1]) * .5f, (min_vert[2] + max_vert[2]) * .5f };

	for (unsigned int i = 0; i < verts.size(); ++i)
	{
		for (unsigned k = 0; k < 3u; ++k)
		{
			verts[i][k] = (verts[i][k] - mid[k]) * fac + mid[k];
		}
	}
}

bool load_obj_verts(std::string const& file_name, Eigen::MatrixXd& V, double scaling_factor,
	std::vector<std::vector<int>> & polygons)
{
	std::vector<std::vector<double>> verts;

	if (!igl::readOBJ(file_name, verts, polygons))
	{
		std::cerr << "Failed to load " << file_name << "!\n";
		return false;
	}

	if (scaling_factor != 1.0)
	{
		scale_verts(verts, scaling_factor);
	}

	V.resize(verts.size(), 3);

	for (int i = 0; i < verts.size(); ++i)
	{
		auto const& vert = verts[i];
		V.row(i) = Eigen::Vector3d(vert[0], vert[1], vert[2]);
	}

	return true;
}

bool load_mesh(std::string const& file_name, Eigen::MatrixXd& V, Eigen::MatrixXi& T, double scaling_factor)
{
	if (file_name.substr(file_name.size() - 4, 4).compare(".msh") == 0)
	{
		if (!igl::readMSH(file_name, V, Eigen::MatrixXi(), T, Eigen::VectorXi(), Eigen::VectorXi()))
		{
			std::cerr << "Failed to load " << file_name << "\n";
			return false;
		}

		if (scaling_factor != 1.0)
		{
			V *= scaling_factor;
		}

	}
	else if (file_name.substr(file_name.size() - 4, 4).compare(".obj") == 0)
	{
		std::vector<std::vector<int>> tris;

		if (!load_obj_verts(file_name, V, scaling_factor, tris))
		{
			return false;
		}

		T.resize(tris.size(), 3);

		for (int i = 0; i < tris.size(); ++i)
		{
			auto const& tri = tris[i];
			assert(tri.size() == 3);
			T.row(i) = Eigen::Vector3i(tri[0], tri[1], tri[2]);
		}

	}
	else
	{
		std::cerr << "Unsupported input format\n";
		return false;
	}

	return true;
}

bool load_cage(std::string const& file_name, Eigen::MatrixXd& V,
	Eigen::VectorXi & P, Eigen::MatrixXi & CF, double scaling_factor, bool triangulate_quads,
	Eigen::MatrixXd * V_embedding /*= nullptr*/, bool find_offset /*= false*/)
{
	if (file_name.substr(file_name.size() - 4, 4).compare(".obj"))
	{
		std::cerr << "Only cages in obj are supported\n";
		return false;
	}

	std::vector<std::vector<int>> polys;

	if (!load_obj_verts(file_name, V, scaling_factor, polys))
	{
		return false;
	}

	int cage_vertices_offset = 0;

	if (V_embedding)
	{
		if (find_offset)
		{
			auto verices_equal = [](const Eigen::Vector3d& a, const Eigen::Vector3d& b)
			{
				return abs(a(0) - b(0)) < igl::FLOAT_EPS && abs(a(1) - b(1)) < igl::FLOAT_EPS && abs(a(2) - b(2)) < igl::FLOAT_EPS;
			};

			const Eigen::Vector3d first_cage_vert = V.row(0);
			bool found = false;
			for (int i = 0; i < V_embedding->rows(); ++i)
			{
				const Eigen::Vector3d embedding_vert = V_embedding->row(i);
				if (verices_equal(first_cage_vert, embedding_vert))
				{
					found = true;
					cage_vertices_offset = i;
					break;
				}
			}
			if (found)
			{
				std::cout << "Found cage verts in embedding with an offset of " << cage_vertices_offset << "\n";
			}
			else
			{
				std::cerr << "Could not find cage verts in embedding\n";
				return 1;
			}
		}
		for (int i = 0; i < V.rows(); ++i)
		{
			V.row(i) = V_embedding->row(i + cage_vertices_offset);
		}
		for (int i = 0; i < CF.rows(); ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				CF(i, j) += cage_vertices_offset;
			}
		}
	}

	P.resize(V.rows());
	P.fill(0);
	for (int i = 0; i < V.rows(); ++i)
	{
		P(i) = i;
	}

	unsigned int numTriangles = 0, numQuads = 0;

	for (int i = 0; i < polys.size(); ++i)
	{
		auto const numVertsPoly = polys[i].size();
		assert(numVertsPoly == 3 || numVertsPoly == 4);

		if (numVertsPoly == 3)
		{
			++numTriangles;
		}
		else // Found a quad
		{
			if (triangulate_quads)
			{
				numTriangles += 2u;
			}
			else
			{
				++numQuads;
			}
		}
	}

	if (triangulate_quads)
	{
		std::cout << "Cage Triangles: " << numTriangles << " Cage Vertices: " << V.rows() << "\n";
	}
	else
	{
		std::cout << "Cage Triangles: " << numTriangles << " Cage Quads " << numQuads << " Cage Vertices: " << V.rows() << "\n";
	}

	CF.resize(triangulate_quads ? numTriangles : numTriangles + numQuads, numQuads ? 4 : 3);

	int row_idx = 0;
	for (int i = 0; i < polys.size(); ++i)
	{
		auto const& poly = polys[i];
		auto const numVertsPoly = poly.size();

		if (numVertsPoly == 3)
		{
			if (numQuads)
			{
				CF.row(row_idx++) = Eigen::Vector4i(poly[0], poly[1], poly[2], -1);
			}
			else
			{
				CF.row(row_idx++) = Eigen::Vector3i(poly[0], poly[1], poly[2]);
			}

		}
		else // quad
		{
			if (triangulate_quads)
			{
				CF.row(row_idx++) = Eigen::Vector3i(poly[0], poly[1], poly[2]);
				CF.row(row_idx++) = Eigen::Vector3i(poly[0], poly[2], poly[3]);
			}
			else
			{
				CF.row(row_idx) = Eigen::Vector4i(poly[0], poly[1], poly[2], poly[3]);
			}
		}
	}

	return true;
}