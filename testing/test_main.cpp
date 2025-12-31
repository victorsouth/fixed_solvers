#define FIXED_USE_QP_SOLVER
#include <fixed/fixed.h>


#define GTEST_BREAK_ON_FAILURE 1
#define GTEST_CATCH_EXCEPTIONS 0
#define GTEST_HAS_SEH 0
#define _VARIADIC_MAX 10 /* for gtest */
#include "gtest/gtest.h"


#include "test_linear_constraints.h"
#include "test_linear_system.h"
#include "test_newton_raphson.h"

#include "qp/eiquadprog.hpp"

Eigen::VectorXd solve_quadprog_box(const Eigen::SparseMatrix<double>& H, const Eigen::VectorXd& f,
    const std::vector<std::pair<size_t, double>>& minimum,
    const std::vector<std::pair<size_t, double>>& maximum
)
{
    Eigen::VectorXd b(minimum.size() + maximum.size());
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(minimum.size() + maximum.size(), f.size());

    size_t row = 0;
    for (const auto& [index, min_bound] : minimum) {
        A(row, index) = -1;
        b(row) = -min_bound;
        row++;
    }
    for (const auto& [index, max_bound] : maximum) {
        A(row, index) = 1;
        b(row) = max_bound;
        row++;
    }


    Eigen::VectorXd estimation = Eigen::VectorXd::Zero(f.size());
    Eigen::MatrixXd CE(estimation.size(), 0); 
    Eigen::VectorXd ce0;
    Eigen::MatrixXd CI = -A.transpose();


    //*result_argmin = initial_estimation;
    double result = solve_quadprog(H, f, CE, ce0, CI, b, estimation);

    return estimation;
}



int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
#ifdef _WIN32
#ifndef __MINGW32__
    std::wcout.imbue(std::locale("rus_rus.866"));
#endif
#endif
int res= RUN_ALL_TESTS();
    return res;
}
