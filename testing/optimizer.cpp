#include "gtest/gtest.h"

#define FIXED_USE_QP_SOLVER
#include <fixed/fixed.h>


/// @brief Проверяет сходимость простой функции
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



/// @brief Проверяет сходимость к верному результату для функции Розенброка
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




auto prepare_ls_task_dense()
{
    // Решаем задачу МНК y = a + b*x
    MatrixXd X(3, 2);
    X << 1, 1,
        1, 2,
        1, 3;

    VectorXd Y(3);
    Y << 2.5, 2, 0.5;

    // Решаем задачу квадратичного программирования, аналогичную МНК
    // Обозначения по матлабу: https://www.mathworks.com/help/optim/ug/quadprog.html
    MatrixXd H = X.transpose() * X;
    VectorXd f = -Y.transpose() * X; // странно, что не надо умножать на 2, но ОК!

    MatrixXd A(1, 2);
    VectorXd b(1);
    A << 1, 0;
    b(0) = 2; // заведомо очень большое ограничение сверху на свободный член


    return std::make_tuple(X, Y, H, f, A, b);
}

SparseMatrix<double, Eigen::ColMajor> dense_to_sparse(const MatrixXd& A) {
    SparseMatrix<double, Eigen::ColMajor> As(A.rows(), A.cols());
    vector<Eigen::Triplet<double>> Atriplets;
    for (int i = 0; i < A.rows(); ++i) {
        for (int j = 0; j < A.cols(); ++j) {
            Atriplets.emplace_back(i, j, A(i, j));
        }
    }
    As.setFromTriplets(Atriplets.begin(), Atriplets.end());
    return As;
}

auto prepare_ls_task_sparse()
{
    // Решаем задачу МНК y = a + b*x
    MatrixXd X(3, 2);
    X << 1, 1,
        1, 2,
        1, 3;

    VectorXd Y(3);
    Y << 2.5, 2, 0.5;

    // Решаем задачу квадратичного программирования, аналогичную МНК
    // Обозначения по матлабу: https://www.mathworks.com/help/optim/ug/quadprog.html
    MatrixXd H = X.transpose() * X;
    VectorXd f = -Y.transpose() * X; // странно, что не надо умножать на 2, но ОК!

    MatrixXd A(1, 2);
    VectorXd b(1);
    A << 1, 0;
    b(0) = 2;

    auto As = dense_to_sparse(A);
    auto Hs = dense_to_sparse(H);

    return std::make_tuple(X, Y, Hs, f, As, b);
}


struct simple_equation_fixed : public fixed_system_t<2> {
    virtual std::array<double, 2> residuals(const std::array<double, 2>& x) {
        std::array<double, 2> x0{ 4, 5 };
        return x - x0;
    }
};


struct simple_equation_var : public fixed_system_t<-1> {
    virtual VectorXd residuals(const VectorXd& x) {
        VectorXd x0(2);
        x0 << 4, 5;
        return x - x0;
    }
};

/// @brief Проверяет способность работы fixed_solvers::newton с квадратичным программированием.
/// Фикс размерность
/// (оставляем здесь)
TEST(Regression, ConstrainedEquations2)
{
    simple_equation_fixed eq;

    fixed_solver_parameters_t<2, 0, golden_section_search> parameters;
    parameters.step_constraint_as_optimization = true;
    parameters.step_constraint_algorithm = step_constraint_algorithm_t::Quadprog;
    //parameters.constraints.minimum = { {0, 7} };
    parameters.constraints.maximum[0] = 3;
    parameters.constraints.minimum[1] = 8;

    std::array<double, 2> x0{ 0, 0 };
    parameters.constraints.ensure_constraints(x0);

    fixed_solver_result_t<2> result;
    fixed_newton_raphson<2>::solve(eq, x0, parameters, &result, nullptr);

    ASSERT_EQ(result.result_code, numerical_result_code_t::Converged);
    ASSERT_NEAR(result.argument[0], parameters.constraints.maximum[0], 1e-8);
    ASSERT_NEAR(result.argument[1], parameters.constraints.minimum[1], 1e-8);

}

/// @brief Проверяет способность работы fixed_solvers::newton с квадратичным программированием.
/// Переменная размерность
/// (оставляем здесь)
TEST(Regression, ConstrainedEquations)
{
    simple_equation_var eq;

    fixed_solver_parameters_t<-1, 0, golden_section_search> parameters;
    parameters.step_constraint_as_optimization = true;
    parameters.step_constraint_algorithm = step_constraint_algorithm_t::Quadprog;
    //parameters.constraints.minimum = { {0, 7} };
    parameters.constraints.maximum = { {0, 3} };

    VectorXd initial = VectorXd::Zero(2);
    parameters.constraints.ensure_constraints(initial);

    fixed_solver_result_t<-1> result;
    fixed_newton_raphson<-1>::solve(eq, initial, parameters, &result);



}

