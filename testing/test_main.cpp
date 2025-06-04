
#define GTEST_BREAK_ON_FAILURE 1
#define GTEST_CATCH_EXCEPTIONS 0
#define GTEST_HAS_SEH 0
#define _VARIADIC_MAX 10 /* for gtest */
#include "gtest/gtest.h"

#define FIXED_USE_QP_SOLVER
#include <fixed/fixed.h>

#include "test_linear_constraints.h"
#include "test_linear_system.h"
#include "test_newton_raphson.h"

#include "qp/eiquadprog.hpp"

VectorXd solve_quadprog_box(const SparseMatrix<double>& H, const VectorXd& f,
    const std::vector<std::pair<size_t, double>>& minimum,
    const std::vector<std::pair<size_t, double>>& maximum
)
{
    VectorXd b(minimum.size() + maximum.size());
    MatrixXd A = MatrixXd::Zero(minimum.size() + maximum.size(), f.size());

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


    VectorXd estimation = VectorXd::Zero(f.size());
    MatrixXd CE(estimation.size(), 0); 
    VectorXd ce0;
    MatrixXd CI = -A.transpose();


    //*result_argmin = initial_estimation;
    double result = solve_quadprog(H, f, CE, ce0, CI, b, estimation);

    return estimation;
}



int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
#ifdef _WIN32
    std::wcout.imbue(std::locale("rus_rus.866"));
#endif
int res= RUN_ALL_TESTS();
    return res;
}
