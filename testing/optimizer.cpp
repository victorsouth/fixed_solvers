#include "gtest/gtest.h"
#include <fixed/fixed.h>



TEST(OptimizeGaussNewton, ConvergesSimpleFunction)
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

TEST(OptimizeGaussNewton, ConvergesRosenbrokFunction)
{
    VectorXd initial = VectorXd::Zero(2);
    rosenbrock_function_t function;

    fixed_solver_parameters_t<-1, 0, golden_section_search> parameters;
    fixed_solver_result_t<-1> result;

    fixed_optimize_gauss_newton::optimize(function, initial, parameters, &result);
    ASSERT_EQ(result.result_code, numerical_result_code_t::Converged);
    ASSERT_NEAR(result.argument(0), 1.0, parameters.argument_increment_norm);
    ASSERT_NEAR(result.argument(1), 1.0, parameters.argument_increment_norm);
};


TEST(OptimizeGaussNewton, PerformsAnalysis)
{
    VectorXd initial = VectorXd::Zero(2);
    rosenbrock_function_t function;

    fixed_solver_parameters_t<-1, 0, golden_section_search> parameters;
    parameters.analysis.objective_function_history = true;
    parameters.analysis.steps = true;
    parameters.analysis.argument_history = true;

    fixed_solver_result_t<-1> result;
    fixed_solver_result_analysis_t<-1> analysis;

    fixed_optimize_gauss_newton::optimize(function, initial, parameters, &result, &analysis);
    ASSERT_EQ(result.result_code, numerical_result_code_t::Converged);

    vector<double> learning_curve = analysis.get_learning_curve();


    FAIL(); // придумать, как проверить
}; 
