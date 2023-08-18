#include <iostream>
#include <string>
#include <array>
#include <vector>
#include <fstream>

#include "GreenCoordinates.h"
#include "InfluenceMap.h"
#include "Parametrization.h"
#include "LoadMesh.h"
#include "WeightInterpolation.h"

#include <boost/program_options.hpp>
//#define VERBOSE
#include <igl/bbw.h>
#include <igl/boundary_conditions.h>
#include <igl/harmonic.h>
#include <igl/lbs_matrix.h>
#include <igl/writeMSH.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/writeDMAT.h>
#include <igl/readDMAT.h>
#include <LBCSolver.h>

int main(int argc, char** argv)
{
	boost::program_options::options_description desc("This is a tool used for morphing with different coordinate variants.\n");

	std::string inputFile, outMeshFile = "a.msh", cageFile, cageDeformedFile, embeddedMeshFile, parameterFile, variant_string = "bbw", modelInfluenceFile;
	int numSamples = 1, numBBWSteps = 20, lbc_scheme = 2;
	float scaling_factor = 1.;

	desc.add_options()
		("model,m", boost::program_options::value<std::string>(&inputFile), "Specifies the input *.btet file of the deformation model")
		("cage,c", boost::program_options::value<std::string>(&cageFile), "Specifies the cage to use (Halfface *.hf file for subspaces and obj for others)")
		("cage-deformed,cd", boost::program_options::value<std::string>(&cageDeformedFile), "Specifies a derformed cage file (instead of a parametrization)")
		("embedded,e", boost::program_options::value<std::string>(&embeddedMeshFile), "Specifies the embedded tet mesh file which embeds the cage (in the first place) and the deformation model (afterwards)")
		("parameters,p", boost::program_options::value<std::string>(&parameterFile), "Specifies the *.param parameter file to control mesh deformation")
		("samples,s", boost::program_options::value<int>(&numSamples), "Specifies the number of samples for the parameters")
		("output,o", boost::program_options::value<std::string>(&outMeshFile), "write resulting mesh to file")
		("no-offset", "Do not calculate an offset for model within the embedding (except the cage)")
		("find-offset", "Search points of embedding for offset (requires equality of points)")
		("scale", boost::program_options::value<float>(&scaling_factor), "Scale model, embedding and cage by factor")
		("iter", boost::program_options::value<int>(&numBBWSteps), "The number of iterations for calculating BBW or LBC")
		("lbc-scheme", boost::program_options::value<int>(&lbc_scheme), "The weighting scheme for lbc")
		("influence", "Evaluate the influence of the control vertices involved in the deformation and write the plot (Mesh for OBJ) to file")
		("interpolate-weights", "Interpolate weights of model vertices from the embedding (embedding does not contain vertices of model)")
		("harmonic", "Use harmonic coordinates by Joshi et al.")
		("LBC", "Use local barycentric coordinates by Zhang et al.")
		("green", "Use green coordinates by Lipman et al.")
		("QGC", "Use tri-quad green coordinates")
		("subspace", "Use Linear subspace design by Wang et al.");
	boost::program_options::positional_options_description p;
	boost::program_options::variables_map vm;
	boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
	boost::program_options::notify(vm);

	const bool harmonic = static_cast<bool>(vm.count("harmonic"));
	const bool lbc = !harmonic && static_cast<bool>(vm.count("LBC"));
	const bool green = !lbc && !harmonic && static_cast<bool>(vm.count("green"));
	const bool QGC = !lbc && !harmonic && !green && static_cast<bool>(vm.count("QGC"));
	const bool find_offset = static_cast<bool>(vm.count("find-offset"));
	const bool scale = static_cast<bool>(vm.count("scale"));
	const bool influence = static_cast<bool>(vm.count("influence"));
	const bool load_deformed_cage = static_cast<bool>(vm.count("cage-deformed"));
	const bool interpolate_weights = static_cast<bool>(vm.count("interpolate-weights"));
	const bool no_offset = static_cast<bool>(vm.count("no-offset")) || interpolate_weights;

	// Display help page if requested
	if (vm.count("help"))
	{
		std::cout << desc << std::endl;
		return 0;
	}

	Parametrization params;
	if (vm.count("parameters"))
	{
		params = readParams(parameterFile);
	}
	else if (!load_deformed_cage)
	{
		std::cerr << "You need to specify either a deformed cage or a parametrization of the cage";
		return 1;
	}

	if (!vm.count("model"))
	{
		std::cerr << "No deformation model specified!\n";
		return 1;
	}

	Eigen::MatrixXd V, V_model;
	Eigen::MatrixXi T, T_model;
	std::cout << "Loading deformation tet mesh\n";
	if (!load_mesh(inputFile, V_model, T_model, scaling_factor))
	{
		return 1;
	}

	if (!green && !QGC && !vm.count("embedded"))
	{
		std::cerr << "You must specify an embedding!\n";
		return 1;
	}

	if (!green && !QGC)
	{
		std::cout << "Loading the embedding\n";
		if (!load_mesh(embeddedMeshFile, V, T, scaling_factor))
		{
			return 1;
		}
	}

	if (!vm.count("cage"))
	{
		std::cerr << "You must specify a cage!\n";
		return 1;
	}

	int model_vertices_offset = 0, cage_vertices_offset = 0;

	// Finding model verts in embedding
	if (!interpolate_weights && find_offset && !green && !QGC)
	{
		auto verices_equal = [](const Eigen::Vector3d& a, const Eigen::Vector3d& b)
		{
			return abs(a(0) - b(0)) < igl::FLOAT_EPS && abs(a(1) - b(1)) < igl::FLOAT_EPS && abs(a(2) - b(2)) < igl::FLOAT_EPS;
		};

		const Eigen::Vector3d first_model_vert = V_model.row(0);
		bool found = false;
		for (int i = 0; i < V.rows(); ++i)
		{
			const Eigen::Vector3d embedding_vert = V.row(i);
			if (verices_equal(first_model_vert, embedding_vert))
			{
				found = true;
				model_vertices_offset = i;
				break;
			}
		}
		if (found)
		{
			std::cout << "Found model verts in embedding with an offset of " << model_vertices_offset << "\n";
		}
		else
		{
			std::cerr << "Could not find model verts in embedding\n";
			return 1;
		}
	}
	else if (!green && !QGC)
	{
		auto const additional_offset = no_offset ? 0 : V.rows() - (V_model.rows() + model_vertices_offset);
		model_vertices_offset += additional_offset;
		std::cout << "Adding an offset of " << additional_offset << "\n";
		std::cout << "Loading embedding\n";
	}

	std::cout << "Loading cage\n";
	Eigen::MatrixXd C, C_deformed;
	Eigen::VectorXi P;
	Eigen::MatrixXi CF, BE, CE;

	if (!load_cage(cageFile, C, P, CF, scaling_factor, !QGC, (!green && !QGC) ? &V : nullptr, find_offset))
	{
		std::cerr << "Failed to load cage!\n";
		return 1;
	}

	if (!find_offset && !interpolate_weights)
	{
		model_vertices_offset = C.rows();
	}

	if (!green && !QGC)
	{
		std::cout << "Using " << model_vertices_offset << " as offset for model vertices in embedding\n";
	}

	if (load_deformed_cage)
	{
		Eigen::VectorXi P_deformed;
		Eigen::MatrixXi CF_deformed;
		if (!load_cage(cageDeformedFile, C_deformed, P_deformed, CF_deformed, scaling_factor, !QGC))
		{
			std::cerr << "Failed to load deformed cage!\n";
			return 1;
		}

		params = fromDeformedCage(C, C_deformed);
	}

	auto const suffix_pos = outMeshFile.find(".");
	const bool write_msh = outMeshFile.substr(suffix_pos + 1, outMeshFile.size()).compare("btet") == 0;

	Eigen::MatrixXd normals;
	std::vector<double> psi_tri;
	std::vector<Eigen::Vector4d> psi_quad;
	Eigen::VectorXi b;
	Eigen::MatrixXd bc;
	if (!green && !QGC)
	{
		std::cout << "Computing Boundary conditions\n";
		if (!igl::boundary_conditions(V, T, C, P, BE, CE, CF, b, bc))
		{
			std::cerr << "Failed to extract boundary conditions for cage!\n";
			return 1;
		}
		std::cout << "Done computing boundary conditions\n";
	}
	else if (green)
	{
		calcNormals(C, CF, normals);
	}
	// compute BBW weights matrix
	Eigen::MatrixXd W, W_interpolated, psi;
	std::cout << "Computing weights\n";
	if (harmonic)
	{
		variant_string = "harmonic";
		if (!igl::harmonic(V, T, b, bc, 1, W))
		{
			std::cerr << "Failed to compute harmonic weights!\n";
			return 1;
		}
	}
	else if (lbc)
	{
		variant_string = std::string("lbc_") + std::to_string(lbc_scheme);
		auto const suffix_pos = embeddedMeshFile.find(".");
		auto const prefix = embeddedMeshFile.substr(0, suffix_pos);
		std::string weightsFile = prefix + std::string("_") + std::to_string(lbc_scheme) + std::string("_lbc_weights.dmat");
		std::ifstream in(weightsFile);
		if (in.good())
		{
			std::cout << "Reading weights from file " << weightsFile << "\n";
			if (!igl::readDMAT(weightsFile, W))
			{
				std::cerr << "Failed to read weights from file!\n";
				return 1;
			}
		}
		else
		{
			LBC::DataSetup::WeightingScheme scheme = static_cast<LBC::DataSetup::WeightingScheme>(lbc_scheme);
			Eigen::MatrixXd sample_points(V.cols(), V.rows());
			for (int i = 0; i < V.rows(); ++i)
			{
				sample_points.col(i) = V.row(i);
			}
			Eigen::MatrixXi cell_vertices(T.cols(), T.rows());
			for (int i = 0; i < T.rows(); ++i)
			{
				cell_vertices.col(i) = T.row(i);
			}
			LBC::IndexVector control_point_idx(P.size());
			for (int i = 0; i < P.size(); ++i)
			{
				control_point_idx(i) = i;
			}
			std::vector< LBC::DataSetup::CageBoundaryFacetInfo > boundary_facet_info;
			for (int i = 0; i < CF.rows(); ++i)
			{
				auto const tri_indices = CF.row(i);
				std::vector<int> boundary_points;

				for (int i = P.size(); i < bc.rows(); ++i)
				{
					auto const row = bc.row(i);
					bool contains = true;
					for (int l = 0; l < 3; ++l)
					{
						if (row(tri_indices(l)) == 0)
						{
							contains = false;
							break;
						}
					}
					if (contains)
					{
						boundary_points.push_back(b(i));
					}
				}

				auto boundary_points_vec = LBC::IndexVector(boundary_points.size());
				for (int i = 0; i < boundary_points.size(); ++i)
				{
					boundary_points_vec(i) = boundary_points[i];
				}

				LBC::DataSetup::CageBoundaryFacetInfo info(tri_indices, boundary_points_vec);
				boundary_facet_info.push_back(info);
			}

			LBC::DataSetup ds(sample_points, control_point_idx, cell_vertices, boundary_facet_info, scheme);

			LBC::Param param;
			param.max_iterations = numBBWSteps;
			param.relaxation_alpha = 1.65;
			param.convergence_check_frequency = 10;
			param.output_frequency_ratio = 10;
			param.rel_primal_eps = 1e-7;
			param.rel_dual_eps = 1e-6;
			param.penalty_weight = 10;
			LBC::LBCSolver solver(param, ds);

			std::cout << "LBC Solver started\n";
			solver.solve();
			std::cout << "Finished computation\n";

			W = ds.get_full_coordinate_values(solver.get_coordinates());
			if (!igl::writeDMAT(weightsFile, W))
			{
				std::cerr << "Failed to write weights to file!\n";
			}
		}
	}
	else if (green)
	{
		variant_string = "green";
		calculateGreenCoordinatesFromQMVC(C, CF, normals, V_model, W, psi);
	}
	else if (QGC)
	{
		variant_string = "QGC";
		calculateGreenCoordinatesTriQuad(C, CF, V_model, W, psi_tri, psi_quad);
	}
	else
	{
		igl::BBWData bbw_data;
		// only a few iterations for sake of demo
		bbw_data.active_set_params.max_iter = numBBWSteps;
		bbw_data.verbosity = 2;
		auto const suffix_pos = embeddedMeshFile.find(".");
		auto const prefix = embeddedMeshFile.substr(0, suffix_pos);
		std::string weightsFile = prefix + std::string("_weights.dmat");
		std::ifstream in(weightsFile);
		if (!in.good())
		{
			if (!igl::bbw(V, T, b, bc, bbw_data, W))
			{
				std::cerr << "Failed to compute bounded biharmonic weights!\n";
				return 1;
			}
			// Write computed weights to file
			if (!igl::writeDMAT(weightsFile, W))
			{
				std::cerr << "Failed to write weights to file!\n";
			}
		}
		else
		{
			std::cout << "Reading weights from file " << weightsFile << "\n";
			if (!igl::readDMAT(weightsFile, W))
			{
				std::cerr << "Failed to read weights from file!\n";
				return 1;
			}
		}
	}
	std::cout << "Done computing weights\n";

	if (!lbc && !green && !QGC)
	{
		igl::normalize_row_sums(W, W);
	}
	Eigen::MatrixXd M;
	// precompute linear blend skinning matrix
	std::cout << "Calculating M\n";
	if (green || QGC)
	{
		M = W;
	}
	else {
		if (interpolate_weights)
		{
			interpolateWeightsInEmbedding(V_model, W, V, T, C.rows(), W_interpolated, M);
		}
		else
		{
			igl::lbs_matrix(V, W, M);
		}
	}
	std::cout << "Done Calculating M\n";
	const int dim = C.cols();
	auto translationFactors = createRegularSampling(numSamples, params.minFactor, params.maxFactor);
	auto const numControlVertices = C.rows();
	std::vector<uint8_t> control_vertices_marked(numControlVertices, 0);
	std::vector<int> control_vertices_idx;
	for (auto it = params.translations_per_vertex.begin(), end = params.translations_per_vertex.end(); it != end; ++it)
	{
		auto const control_vertex_idx = it->first;

		control_vertices_marked[control_vertex_idx] = 1;
		if (influence)
		{
			control_vertices_idx.push_back(control_vertex_idx);
		}
	}

	if (influence)
	{
		auto const base_name = outMeshFile.substr(0, suffix_pos);
		write_influence_color_map_OBJ(base_name + "_influence_" + variant_string + ".obj", V_model, T_model, interpolate_weights ? W_interpolated : W,
			control_vertices_idx, (green || QGC || interpolate_weights) ? 0 : model_vertices_offset, green || QGC);
	}

	auto const numTransformations = numControlVertices;
	Eigen::VectorXi tet_tags(T_model.rows());
	for (int j = 0; j < T_model.rows(); ++j)
	{
		tet_tags(j) = 1;
	}

	Eigen::MatrixXd Transformation(numTransformations * (dim + 1), dim);
	for (int j = 0; j < numTransformations; ++j)
	{
		auto a = Eigen::Affine3d::Identity();
		Transformation.block(j * (dim + 1), 0, dim + 1, dim) =
			a.matrix().transpose().block(0, 0, dim + 1, dim);
	}

	for (int i = 0; i < translationFactors.size(); ++i)
	{
		Eigen::MatrixXd Cage_transforms = Transformation;
		C_deformed = C;
		auto const factor = translationFactors[i];

		for (auto it = params.translations_per_vertex.begin(), end = params.translations_per_vertex.end(); it != end; ++it)
		{
			auto const cage_vertex_idx = it->first;
			auto const translation = Eigen::Affine3d(Eigen::Translation3d(factor * it->second));
			Transformation.block(cage_vertex_idx * (dim + 1), 0, dim + 1, dim) =
				translation.matrix().transpose().block(0, 0, dim + 1, dim);
			Eigen::Vector3d cage_vertex = C.row(cage_vertex_idx);
			C_deformed.row(cage_vertex_idx) = cage_vertex + factor * it->second;
		}

		Eigen::MatrixXd U, U_model(V_model.rows(), 3);
		if (green)
		{
			Eigen::MatrixXd normals_deformed;
			calcNormals(C_deformed, CF, normals_deformed);
			calcScalingFactors(C, C_deformed, CF, normals_deformed);
			U_model = W.transpose() * C_deformed + psi.transpose() * normals_deformed;
		}
		else if (QGC)
		{
			calcNewPositionsTriQuad(C, C_deformed, CF, W, psi_tri, psi_quad, U_model);
		}
		else
		{
			U = M * Transformation;
		}

		if (!green && !QGC)
		{
			for (int j = 0; j < V_model.rows(); ++j)
			{
				U_model.row(j) = U.row(j + model_vertices_offset);
				if (i == 0)
				{
					Eigen::Vector3d a = U_model.row(j);
					auto const vert = V_model.row(j);
					Eigen::Vector3d b(vert[0], vert[1], vert[2]);
					const float dist = (a - b).squaredNorm();
					assert(dist < igl::FLOAT_EPS);
				}

			}
		}

		auto const prefix = outMeshFile.substr(0, suffix_pos);
		auto const interpolation_factor = translationFactors.size() > 1 ? static_cast<float>(i) / static_cast<float>(translationFactors.size() - 1) : 1;
		auto const middle = std::string("_") + std::to_string(interpolation_factor) + std::string("_");

		//igl::writeOBJ(prefix + std::string("_cage") + middle + variant_string + std::string(".obj"), C_deformed, CF);
		if (write_msh)
		{
			igl::writeMSH(prefix + middle + variant_string + std::string(".msh"), U_model, Eigen::MatrixXi(), T_model, Eigen::MatrixXi(), tet_tags, std::vector<std::string>(),
				std::vector<Eigen::MatrixXd>(), std::vector<std::string>(), std::vector<Eigen::MatrixXd>(), std::vector<Eigen::MatrixXd>());
		}
		else
		{
			std::cout << "Writing " << prefix + middle + variant_string + std::string(".obj") << "\n";
			igl::writeOBJ(prefix + middle + variant_string + std::string(".obj"), U_model, T_model);
		}
	}

	return 0;
}
