#include "gtest/gtest.h"

#define FIXED_USE_QP_SOLVER
#include <fixed/fixed.h>

/// @brief Простейшая система уравнений 2x2 фиксированной размерности
struct simple_equation_fixed : public fixed_system_t<2> {
    virtual std::array<double, 2> residuals(const std::array<double, 2>& x) {
        std::array<double, 2> x0{ 4, 5 };
        return x - x0;
    }
};

/// @brief Простейшая система уравнений 2x2 размерность переменная
struct simple_equation_var : public fixed_system_t<-1> {
    virtual VectorXd residuals(const VectorXd& x) {
        VectorXd x0(2);
        x0 << 4, 5;
        return x - x0;
    }
};

/// @brief Проверяет способность работы fixed_solvers::newton с квадратичным программированием.
/// Фиксированная размерность
TEST(NewtonRaphson, HandlesConstrainedEquationsFixed)
{
    simple_equation_fixed eq;

    fixed_solver_parameters_t<2, 0, golden_section_search> parameters;
    parameters.step_constraint_as_optimization = true;
    parameters.step_constraint_algorithm = step_constraint_algorithm_t::Quadprog;
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
TEST(NewtonRaphson, HandlesConstrainedEquationsVar)
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

    ASSERT_EQ(result.result_code, numerical_result_code_t::Converged);
    ASSERT_NEAR(result.argument(0), parameters.constraints.maximum.front().second, 1e-8);
    ASSERT_NEAR(result.residuals(1), 0.0, 1e-8);

}

