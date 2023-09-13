#define _USE_MATH_DEFINES
#include <cmath>

#include "MaximumEntropyCoordinates.h"

#include <iostream>
#include <vector>
#include <Eigen/Geometry>

void calculateMaximumEntropyCoordinates(const Eigen::MatrixXd &cage_v, const Eigen::MatrixXi &cage_f, const Eigen::MatrixXd &model_v, Eigen::MatrixXd &mec, const int mec_flag)
{
    int nv = cage_v.cols();
    int nf = cage_f.cols();

    Eigen::MatrixXd shifted_cage_v = cage_v;

    mec.resize(nv, model_v.cols());

    std::vector<std::vector<int>> allAdjFaces;
    std::vector<int> adjFaces;
    // get adjacent faces of each vertex
    for (int i = 0; i < nv; ++i)
    {
        adjFaces = findAdjacentFacesOfIndexV(cage_f, i);
        allAdjFaces.push_back(adjFaces);
    }

    for (int ii = 0; ii < model_v.cols(); ++ii)
    {

        // compute the prior function
        Eigen::VectorXd priors = Eigen::VectorXd::Ones(nv);
        priorFunctions(model_v.col(ii), cage_v, cage_f, allAdjFaces, priors, mec_flag);

        for (int i = 0; i < nv; ++i)
        {
            shifted_cage_v.col(i) = cage_v.col(i) - model_v.col(ii);
        }

        // damp Newton
        double error = 1e-10;
        double rho = .55;
        double sigma = .4;
        int mm, mk;
        Eigen::Vector3d lambda(0.0, 0.0, 0.0);
        while (true)
        {
            Eigen::VectorXd Zi = Eigen::VectorXd::Zero(nv);
            for (int i = 0; i < nv; ++i)
            {
                Zi[i] = priors[i] * exp(-lambda.dot(shifted_cage_v.col(i)));
            }
            // compute the gradient
            Eigen::Vector3d gradient(0.0, 0.0, 0.0);
            for (int i = 0; i < nv; ++i)
            {
                gradient -= Zi[i] * shifted_cage_v.col(i);
            }
            double Z = Zi.sum();
            gradient /= Z;

            if (gradient.norm() < error)
            {
                mec.col(ii) = Zi / Z;
                // std::cout << (cage_v * mec.col(ii) - model_v.col(ii)).norm() << std::endl; // for test
                break;
            }

            // compute the Hessian matrix
            Eigen::Vector3d SUMZivi(0.0, 0.0, 0.0);
            Eigen::Matrix3d Hessian = Eigen::Matrix3d::Zero();
            for (int i = 0; i < nv; ++i)
            {
                Eigen::Vector3d Zivi = Zi[i] * shifted_cage_v.col(i);
                SUMZivi += Zivi;
                Hessian += Zivi * (shifted_cage_v.col(i).transpose()); // 3-by-3
            }
            Hessian *= Z;
            Hessian -= SUMZivi * SUMZivi.transpose();
            Hessian /= (Z * Z);

            // determine Newton search direction
            SUMZivi = -Hessian.inverse() * gradient;

            // determine step size by using Armijo linear search
            mm = 0;
            mk = 0;
            while (mm < 20)
            {
                Eigen::Vector3d lambda_temp = lambda + pow(rho, mm) * SUMZivi;
                for (int i = 0; i < nv; ++i)
                {
                    Zi[i] = priors[i] * exp(-lambda_temp.dot(shifted_cage_v.col(i)));
                }
                if (log(Zi.sum()) < log(Z) + sigma * pow(rho, mm) * gradient.dot(SUMZivi))
                {
                    mk = mm;
                    break;
                }
                mm++;
            }
            lambda += pow(rho, mk) * SUMZivi; // update lambda
        }
    }
}

std::vector<int> findAdjacentFacesOfIndexV(const Eigen::MatrixXi &faces, int indexOfVertex)
{
    std::vector<int> adjacentFaces;
    for (int i = 0; i < faces.cols(); ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            if (faces(j, i) == indexOfVertex)
            {
                adjacentFaces.push_back(i);
                break;
            }
        }
    }
    return adjacentFaces;
}

void priorFunctions(const Eigen::Vector3d v, const Eigen::MatrixXd &cage_v, const Eigen::MatrixXi &cage_f, std::vector<std::vector<int>> adjs, Eigen::VectorXd &priors, int mec_flag)
{
    for (int i = 0; i < cage_v.cols(); ++i)
    {

        if (mec_flag == 1) // MEC-1
        {
            for (std::vector<int>::iterator iter = adjs.at(i).begin(); iter != adjs.at(i).end(); ++iter)
            {
                priors[i] *= areaOfTriangle(v, cage_v.col(cage_f(0, *iter)), cage_v.col(cage_f(1, *iter))) + areaOfTriangle(v, cage_v.col(cage_f(1, *iter)), cage_v.col(cage_f(2, *iter))) +
                             areaOfTriangle(v, cage_v.col(cage_f(2, *iter)), cage_v.col(cage_f(0, *iter))) - areaOfTriangle(cage_v.col(cage_f(0, *iter)), cage_v.col(cage_f(1, *iter)), cage_v.col(cage_f(2, *iter)));
            }
        }
        if (mec_flag == 2) // MEC-2
        {
            for (std::vector<int>::iterator iter = adjs.at(i).begin(); iter != adjs.at(i).end(); ++iter)
            {
                priors[i] *= areaOfTriangle(v, cage_v.col(cage_f(0, *iter)), cage_v.col(cage_f(1, *iter))) * areaOfTriangle(v, cage_v.col(cage_f(1, *iter)), cage_v.col(cage_f(2, *iter))) *
                                 areaOfTriangle(v, cage_v.col(cage_f(2, *iter)), cage_v.col(cage_f(0, *iter))) +
                             (v - cage_v.col(cage_f(0, *iter))).dot(v - cage_v.col(cage_f(1, *iter))) + (v - cage_v.col(cage_f(1, *iter))).dot(v - cage_v.col(cage_f(2, *iter))) + (v - cage_v.col(cage_f(2, *iter))).dot(v - cage_v.col(cage_f(0, *iter)));
            }
        }
        priors[i] = 1.0 / priors[i];
    }
    priors /= priors.sum();
}

double areaOfTriangle(const Eigen::Vector3d a, const Eigen::Vector3d b, const Eigen::Vector3d c)
{
    return 0.5 * (b - a).cross(c - a).norm();
}