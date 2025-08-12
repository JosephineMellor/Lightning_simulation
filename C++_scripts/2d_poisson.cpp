#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <array>
#include <filesystem>
#include <tuple>
#include <Eigen/Sparse>

//set up the grid to use for the finite difference scheme
const int Nr = 100; 
const int Nz = 100;
const int N = Nr * Nz;

//set up ranges of r and z and find differences
const double r0 = 0.0, r1 = 1.0;
const double z0 = 0.0, z1 = 1.0;
const int I0 = 106405;
const int alpha = 22708;
const int beta = 1294530;
const int gamma = 10847100;
const double PI = 3.141592653589793;
const double r0 = 2e-3; //2mm in meters

//flatten the (i,j) into the index of a vector of size N (there will be Nz of each i with each j added on)
int vectoridx(int i, int j){
    return i * Nz + j;
}

//set up current distribution
std::array<double,2> current(double r, double z, double t){
    double I = I0 * (std::exp(-alpha * t) - std::exp(-beta * t)) * ( 1 - std::exp(-gamma * t)) * ( 1 - std::exp(-gamma * t));
    std::array<double,2> J;
    J[1] = - (I/PI*r0*r0) * std::exp(-(r/r0)*(r/r0));
    J[0] = 0;
    return J;
}

//solve the sparse matrix linear system
std::array<Eigen::VectorXd,2> poissonSolver(double t){
    double dr = (r1 - r0) / (Nr - 1);
    double dz = (z1 - z0) / (Nz - 1);
    double dr2 = dr * dr;
    double dz2 = dz * dz;

    Eigen::SparseMatrix<double> M(N,N);
    Eigen::VectorXd bR = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd bZ = Eigen::VectorXd::Zero(N);
    std::vector<Eigen::Triplet<double>> coefficients;

    //set up matrix
    for(int i =0; i < Nr; ++i){
        double r = i * dr;

        for(int j =0; j < Nz; ++j){
            int k = vectoridx(i,j);
            double z = j * dz;

            //dirichlet conditions
            if(i == 0 || i == Nr - 1 || j == 0 || j == Nz - 1){
                coefficients.emplace_back(k,k,1.0);
                bR[k] = 0.0;
                bZ[k] = 0.0;
                continue;
            }

            //coefficients for finite difference scheme for poisson equation
            double rPlus = r + dr / 2.0;
            double rMinus = r - dr / 2.0;
            double centre = -2.0 / dr2 - 2.0 / dz2;
            double coeffPlusr = (1.0 + dr / (2.0*r)) / dr2;
            double coeffMinusr = (1.0 - dr / (2.0*r)) / dr2;
            double coeffz = 1.0 / dz2;

            // Fill matrix entries
            coefficients.emplace_back(k, vectoridx(i, j), centre);
            coefficients.emplace_back(k, vectoridx(i + 1, j), coeffPlusr);
            coefficients.emplace_back(k, vectoridx(i - 1, j), coeffMinusr);
            coefficients.emplace_back(k, vectoridx(i, j + 1), coeffz);
            coefficients.emplace_back(k, vectoridx(i, j - 1), coeffz);

            //RHS
            bR[k] = -current(r,z,t)[0];
            bZ[k] = -current(r,z,t)[1];
        }
    }

    //build a sparse matrix
    M.setFromTriplets(coefficients.begin(), coefficients.end());

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(M);
    solver.factorize(M);
    Eigen::VectorXd AR = solver.solve(bR);
    Eigen::VectorXd AZ = solver.solve(bZ);

    std::array<Eigen::VectorXd,2> A;
    A[0] = AR;
    A[1] = AZ;

    return A;
}

