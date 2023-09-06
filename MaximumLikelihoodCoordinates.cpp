#define _USE_MATH_DEFINES
#include <cmath>

#include "MaximumLikelihoodCoordinates.h"

#include <iostream>
#include <vector>
#include <Eigen/Geometry>

void computeIntegralUnitNormals(const Eigen::MatrixXd &v, const Eigen::MatrixXi &f, Eigen::MatrixXd &transMatrix, Eigen::MatrixXd &integral_outward_allfaces)
{
    integral_outward_allfaces.resize(3, f.cols());
    Eigen::Matrix3d T, norm_;
    Eigen::Vector3d cross_0, cross_1, cross_2, length_c, beta;
    for (int i = 0; i < f.cols(); ++i)
    {
        T.col(0) = v.col(f(0, i));
        T.col(1) = v.col(f(1, i));
        T.col(2) = v.col(f(2, i));

        norm_.col(0) = T.col(0).cross(T.col(1));
        norm_.col(0) = norm_.col(0).normalized();
        norm_.col(1) = T.col(1).cross(T.col(2));
        norm_.col(1) = norm_.col(1).normalized();
        norm_.col(2) = T.col(2).cross(T.col(0));
        norm_.col(2) = norm_.col(2).normalized();

        beta[0] = acos(T.col(0).dot(T.col(1)));
        beta[1] = acos(T.col(1).dot(T.col(2)));
        beta[2] = acos(T.col(2).dot(T.col(0)));

        integral_outward_allfaces.col(i) = norm_ * beta / 2.0;

        for (int j = 0; j < 3; ++j)
        {
            transMatrix(j, i) =
                (beta[(j + 1) % 3] + beta[j] * norm_.col(j).dot(norm_.col((j + 1) % 3)) +
                 beta[(j + 2) % 3] * norm_.col((j + 2) % 3).dot(norm_.col((j + 1) % 3))) /
                (2 * T.col(j).dot(norm_.col((j + 1) % 3)));
        }
        if (T.determinant() < 0)
        {
            integral_outward_allfaces.col(i) *= -1;
            transMatrix.col(i) *= -1;
        }
    }
}

void calculateMaximumLikelihoodCoordinates(const Eigen::MatrixXd &cage_v, const Eigen::MatrixXi &cage_f, const Eigen::MatrixXd &model_v, Eigen::MatrixXd &mlc)
{
    int nv = cage_v.cols();
    int nf = cage_f.cols();

    mlc.resize(nv, model_v.cols());

    std::vector<std::vector<int>> allAdjFaces;
    std::vector<int> adjFaces;
    // get adjacent faces of each vertex
    for (int i = 0; i < nv; ++i)
    {
        adjFaces = findAdjacentFaces(cage_f, i);
        allAdjFaces.push_back(adjFaces);
    }

    for (int ii = 0; ii < model_v.cols(); ++ii)
    {
        // Eigen::VectorXd mlc = Eigen::VectorXd::Zero(nv);
        std::vector<Eigen::MatrixXd> transMatrixGroup;
        Eigen::MatrixXd transMatrix = Eigen::MatrixXd::Identity(nv, nv);

        // 1. get projection vertices on unit sphere
        Eigen::MatrixXd v;
        v.resize(3, nv);
        for (int i = 0; i < nv; ++i)
        {
            v.col(i) = cage_v.col(i) - model_v.col(ii);
            double r = v.col(i).norm();
            transMatrix(i, i) = 1.0 / r;
        }
        transMatrixGroup.push_back(transMatrix);
        v *= transMatrix;
        // std::cout << v.transpose() << std::endl;

        // 2. smooth
        Eigen::MatrixXd tempVar;
        double integral_vector_length;

        tempVar = Eigen::MatrixXd::Zero(3, nf);
        transMatrix = Eigen::MatrixXd::Zero(nv, nv);

        Eigen::MatrixXd integral_outward_allfaces;
        computeIntegralUnitNormals(v, cage_f, tempVar, integral_outward_allfaces);

        v = Eigen::MatrixXd::Zero(3, nv);

        for (int i = 0; i < nv; ++i)
        {
            for (std::vector<int>::iterator iter = allAdjFaces.at(i).begin();
                 iter != allAdjFaces.at(i).end(); ++iter)
            {
                v.col(i) += integral_outward_allfaces.col(*iter);

                transMatrix(cage_f(0, *iter), i) += tempVar(0, *iter);
                transMatrix(cage_f(1, *iter), i) += tempVar(1, *iter);
                transMatrix(cage_f(2, *iter), i) += tempVar(2, *iter);
            }
            integral_vector_length = v.col(i).norm();
            transMatrix.col(i) /= integral_vector_length;
            v.col(i) /= integral_vector_length;
        }

        transMatrixGroup.push_back(transMatrix);

        // damp Newton
        Eigen::Vector3d x(0, 0, 0);
        Eigen::VectorXd g, d;
        Eigen::Matrix3d H;
        double error = 1e-10;
        double rho = .55;
        double sigma = .4;
        int mm, mk;
        while (true)
        {
            g = f_gradient(x, v);
            if (g.norm() < error)
            {
                break;
            }
            H = f_Hessian(x, v);
            d = -H.inverse() * g;

            // Armijo linear search
            mm = 0;
            mk = 0;
            while (mm < 20)
            {
                if (f(x + pow(rho, mm) * d, v) < f(x, v) + sigma * pow(rho, mm) * g.dot(d))
                {
                    mk = mm;
                    break;
                }
                mm++;
            }
            // cout << pow(rho, mk) << endl;
            x += pow(rho, mk) * d;
        }

        for (int i = 0; i < nv; ++i)
        {
            mlc(i, ii) = 1 / (nv + x[0] * v(0, i) + x[1] * v(1, i) + x[2] * v(2, i));
        }

        // std::cout << mlc << std::endl;
        while (!transMatrixGroup.empty())
        {
            mlc.col(ii) = transMatrixGroup.back() * mlc.col(ii);
            transMatrixGroup.pop_back();
        }

        mlc.col(ii) /= mlc.col(ii).sum();
    }
}

std::vector<int> findAdjacentFaces(const Eigen::MatrixXi &faces, int indexOfVertex)
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

Eigen::VectorXd f_gradient(Eigen::VectorXd x, Eigen::MatrixXd v)
{
    int n = v.cols();
    Eigen::Vector3d g(0, 0, 0);

    for (int i = 0; i < n; ++i)
    {
        g -= v.col(i) / (n + x.dot(v.col(i)));
    }

    return g;
}

Eigen::MatrixXd f_Hessian(Eigen::VectorXd x, Eigen::MatrixXd v)
{
    int n = v.cols();

    Eigen::Matrix3d H = Eigen::Matrix3d::Zero();

    for (int i = 0; i < 3; ++i)
    {
        for (int j = i; j < 3; ++j)
        {
            for (int k = 0; k < n; ++k)
            {
                H(i, j) += v(i, k) * v(j, k) / ((n + x.dot(v.col(k))) * (n + x.dot(v.col(k))));
            }
            H(j, i) = H(i, j);
        }
    }
    return H;
}

double f(Eigen::VectorXd x, Eigen::MatrixXd v)
{
    int n = v.cols();
    double z = 0;

    for (int i = 0; i < n; ++i)
    {
        z -= log(n + x.dot(v.col(i)));
    }

    return z;
}