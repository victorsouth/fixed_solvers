#define FIXED_USE_QP_SOLVER
#include <fixed/fixed.h>

#include "gtest/gtest.h"


/// @brief Простейшая система уравнений 2x2 фиксированной размерности
struct simple_equation_fixed : public fixed_system_t<2> {
    /// @brief Невязки
    virtual std::array<double, 2> residuals(const std::array<double, 2>& x) {
        std::array<double, 2> x0{ 4, 5 };
        return x - x0;
    }
};

/// @brief Простейшая система уравнений 2x2 размерность переменная
struct simple_equation_var : public fixed_system_t<-1> {
    /// @brief Невязки
    virtual Eigen::VectorXd residuals(const Eigen::VectorXd& x) {
        Eigen::VectorXd x0(2);
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
/// Фиксированная размерность
TEST(NewtonRaphson, HandlesConstrainedEquationsCoordinateDescent)
{
    simple_equation_fixed eq;

    fixed_solver_parameters_t<2, 0, golden_section_search> parameters;
    parameters.step_constraint_as_optimization = true;
    parameters.step_constraint_algorithm = step_constraint_algorithm_t::CoordinateDescent;
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

    Eigen::VectorXd initial = Eigen::VectorXd::Zero(2);
    parameters.constraints.ensure_constraints(initial);

    fixed_solver_result_t<-1> result;
    fixed_newton_raphson<-1>::solve(eq, initial, parameters, &result);

    ASSERT_EQ(result.result_code, numerical_result_code_t::Converged);
    ASSERT_NEAR(result.argument(0), parameters.constraints.maximum.front().second, 1e-8);
    ASSERT_NEAR(result.residuals(1), 0.0, 1e-8);

}

/// @brief Простая система переменной размерности с диагностикой исследования шага
struct diagnostic_equation_var : public fixed_system_t<-1> {
    /// @brief Собранные значения alpha по шагам Ньютона-Рафсона
    std::vector<std::vector<double>> alphas_per_step;

    /// @brief Невязки
    virtual Eigen::VectorXd residuals(const Eigen::VectorXd& x) override {
        Eigen::VectorXd x0(2);
        x0 << 1.0, -2.0;
        return x - x0;
    }

    /// @brief Пользовательское исследование траектории шага - запоминает сетку alpha
    virtual void custom_line_research(const Eigen::VectorXd& argument, const Eigen::VectorXd& argument_increment) override {
        (void)argument;
        (void)argument_increment;
        alphas_per_step.emplace_back();
        auto& alphas = alphas_per_step.back();
        size_t research_step_count = 100;
        for (size_t index = 0; index <= research_step_count; ++index) {
            double alpha = 1.0 * index / research_step_count;
            alphas.push_back(alpha);
        }
    }
};

/// @brief Проверяет возможность сбора базовой и пользовательской диагностики 
/// при исследовании траектории шага метода Ньютона-Рафсона
TEST(NewtonRaphson, HandlesLineSearchCustomDiagnostics)
{
    // Arrange - Создание системы с пользовательской диагностикой. Включение диагностики.
    diagnostic_equation_var eq;

    fixed_solver_parameters_t<-1, 0, golden_section_search> parameters;
    parameters.analysis.line_search_explore = true;

    Eigen::VectorXd initial = Eigen::VectorXd::Zero(2);

    fixed_solver_result_t<-1> result;
    fixed_solver_result_analysis_t<-1> analysis;

    // Act - Расчет с активным флагом исследования
    fixed_newton_raphson<-1>::solve(eq, initial, parameters, &result, &analysis);

    // Arrange - должен быть хотя бы один шаг с исследованием
    ASSERT_FALSE(analysis.target_function.empty());
    ASSERT_FALSE(eq.alphas_per_step.empty());

    // Количество шагов совпадает в общей и в кастомной диагностике
    ASSERT_EQ(analysis.target_function.size(), eq.alphas_per_step.size());

    // Для первого шага количество элементов в общей и кастомной диагностике совпадает
    ASSERT_EQ(analysis.target_function.front().size(), eq.alphas_per_step.front().size());
    ASSERT_GT(eq.alphas_per_step.front().size(), 0u);
}

/// @brief Простое кубическое уравнение
struct cubic_equation_fixed : public fixed_system_t<1> {
    /// @brief Невязки
    virtual double residuals(const double& x) override {
        double x0 = 3;
        double r = std::pow(x - x0, 3.0);
        return r;
    }
    /// @brief Переопределяем целевую функцию, чтобы был модуль невязок
    virtual double objective_function(const double& r) const override {
        return std::abs(r);
    }
};

/// @brief Проверят способность учесть ограничения по точности по невязке
TEST(NewtonRaphson, HandlesResidualsNorm) {

    // Arrange
    cubic_equation_fixed equation;

    fixed_solver_parameters_t<1, 0, golden_section_search> parameters;
    parameters.residuals_norm = 1.5; // допускаем невязку целевой функции до 1.5
    parameters.residuals_norm_allow_early_exit = true;

    // Act
    // начальное приближение берем далеко от решения x = 3.0, чтобы изначально была невязка
    double initial_x = 10; 
    fixed_solver_result_t<1> result;
    fixed_newton_raphson<1>::solve(equation, initial_x, parameters, &result);

    // Assert - проверяем, что процесс сошелся, но невязкам, а не по приращению аргумента
    double norm = equation(result.argument);
    ASSERT_LE(norm, parameters.residuals_norm); // невязки действительно не превышают
    ASSERT_EQ(result.result_code, numerical_result_code_t::Converged); // сошелся
    ASSERT_FALSE(result.argument_increment_criteria); // критерий приращения не удовлетоврен
    ASSERT_TRUE(result.residuals_norm_criteria); // критерий по невязке удовлетворен
    

}
