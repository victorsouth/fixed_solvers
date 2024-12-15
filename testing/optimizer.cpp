#include "gtest/gtest.h"
#include <fixed/fixed.h>


/// @brief �������� ����������� ��������� ������� �������
/// J(x1, x2) = (x1 - 2)^2 + (x2 - 1)^2 
TEST(OptimizeGaussNewton, ConvergesSimpleFunction)
{
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

    fixed_optimizer_parameters_t parameters;
    fixed_optimizer_result_t result;

    fixed_optimize_gauss_newton::optimize(function, initial, parameters, &result);
    ASSERT_EQ(result.result_code, numerical_result_code_t::Converged);
    ASSERT_NEAR(result.argument(0), 2.0, parameters.argument_increment_norm);
    ASSERT_NEAR(result.argument(1), 1.0, parameters.argument_increment_norm);
}



/// @brief ��������� ����������� �������������� ������� ����������
TEST(OptimizeGaussNewton, ConvergesRosenbrokFunction)
{
    VectorXd initial = VectorXd::Zero(2);
    rosenbrock_function_t function;

    fixed_optimizer_parameters_t parameters;
    fixed_optimizer_result_t result;

    fixed_optimize_gauss_newton::optimize(function, initial, parameters, &result);
    ASSERT_EQ(result.result_code, numerical_result_code_t::Converged);
    ASSERT_NEAR(result.argument(0), 1.0, parameters.argument_increment_norm);
    ASSERT_NEAR(result.argument(1), 1.0, parameters.argument_increment_norm);
};

/// @brief Проверяет способность собирать кривую обучение 
/// Проверяется снижение целевой функции на каждом шаге расчета
TEST(OptimizeGaussNewton, PerformsLearningCurveAnalysis)
{
    VectorXd initial = VectorXd::Zero(2);
    rosenbrock_function_t function;

    fixed_optimizer_parameters_t parameters;
    parameters.analysis.objective_function_history = true;
    parameters.analysis.steps = true;
    parameters.analysis.argument_history = true;

    fixed_optimizer_result_t result;
    fixed_optimizer_result_analysis_t analysis;

    fixed_optimize_gauss_newton::optimize(function, initial, parameters, &result, &analysis);
    ASSERT_EQ(result.result_code, numerical_result_code_t::Converged);

    vector<double> learning_curve = analysis.get_learning_curve();

    ASSERT_EQ(learning_curve.empty(), false);

    for (size_t index = 1; index < learning_curve.size(); ++index) {
        bool decrements = learning_curve[index] < learning_curve[index - 1];
        ASSERT_EQ(decrements, true);
    }
}; 
