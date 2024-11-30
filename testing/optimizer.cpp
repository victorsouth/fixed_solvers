#include "gtest/gtest.h"
#include <fixed/fixed.h>



TEST(OptimizeGaussNewton, Test1)
{
    /// @brief J(x1, x2) = (x1 - 2)^2 + (x2 - 1)^2 
    class simple_sum_of_squares_function : public fixed_least_squares_function_t
    {
    public:
        VectorXd residuals(const VectorXd& x) {
            VectorXd result(2);
            result[0] = x(0) - 2.0;
            result[1] = x(1) - 1.0;
            return result;
        }
    };

    VectorXd initial = VectorXd::Zero(2);
    simple_sum_of_squares_function function;

    fixed_solver_parameters_t<-1, 0, golden_section_search> parameters;
    fixed_solver_result_t<-1> result;

    fixed_optimize_gauss_newton::optimize(function, initial, parameters, &result);
    ASSERT_EQ(result.result_code, numerical_result_code_t::Converged);
    ASSERT_NEAR(result.argument(0), 2.0, parameters.argument_increment_norm);
    ASSERT_NEAR(result.argument(1), 1.0, parameters.argument_increment_norm);
}

TEST(OptimizeGaussNewton, RosenbrokFunction)
{
    /// @brief J(x1, x2) = (x1 - 2)^2 + (x2 - 1)^2 
    class rosenbrock_function_t : public fixed_least_squares_function_t
    {
    public:
        VectorXd residuals(const VectorXd& x) {
            VectorXd result(2);
            result[0] = 10 * (x(1) - x(0) * x(0));
            result[1] = 1 - x(0);
            return result;
        }
    };

    VectorXd initial = VectorXd::Zero(2);
    rosenbrock_function_t function;

    fixed_solver_parameters_t<-1, 0, golden_section_search> parameters;
    fixed_solver_result_t<-1> result;

    fixed_optimize_gauss_newton::optimize(function, initial, parameters, &result);
    ASSERT_EQ(result.result_code, numerical_result_code_t::Converged);
    ASSERT_NEAR(result.argument(0), 1.0, parameters.argument_increment_norm);
    ASSERT_NEAR(result.argument(1), 1.0, parameters.argument_increment_norm);
};
