#pragma once

namespace {
/// @brief Система: ООФ x>=0 и при конечном upper — x<=upper
struct residual_domain_violation_var : public fixed_system_t<-1> {
    double upper_bound_for_domain{ std::numeric_limits<double>::infinity() };
    double linear_residual_root{ 1.0 };
    Eigen::VectorXd residuals(const Eigen::VectorXd& x) override {
        if (x(0) < 0.0) {
            throw domain_violation{};
        }
        if (std::isfinite(upper_bound_for_domain) && x(0) > upper_bound_for_domain) {
            throw domain_violation{};
        }
        Eigen::VectorXd r(1);
        r(0) = x(0) - linear_residual_root;
        return r;
    }
};

/// @brief Система переменной размерности, кидающая domain_violation при расчете Якобиана.
struct jacobian_domain_violation_var : public fixed_system_t<-1> {
    Eigen::VectorXd residuals(const Eigen::VectorXd& x) override {
        Eigen::VectorXd r(1);
        r(0) = x(0) - 1.0;
        return r;
    }
    sparse_matrix_type jacobian_sparse(const Eigen::VectorXd& x) override {
        throw domain_violation{};
    }
};

/// @brief Система фиксированной размерности 2, кидающая domain_violation в residuals.
struct residual_domain_violation_fixed2 : public fixed_system_t<2> {
    std::array<double, 2> residuals(const std::array<double, 2>& x) override {
        if (x[0] < 0.0) {
            throw domain_violation{};
        }
        return { x[0] - 1.0, x[1] - 2.0 };
    }
};

/// @brief Система фиксированной размерности 2, кидающая domain_violation в jacobian_dense.
struct jacobian_domain_violation_fixed2 : public fixed_system_t<2> {
    std::array<double, 2> residuals(const std::array<double, 2>& x) override {
        return { x[0] - 1.0, x[1] - 2.0 };
    }
    matrix_value jacobian_dense(const var_type& x) override {
        throw domain_violation{};
    }
};

}

/// @brief Проверяет перехват domain_violation при недопустимом начальном аргументе для fixed_system_t<-1>.
TEST(NewtonRaphsonDomainViolation, CatchesResidualViolationAtInitialPointVar)
{
    // Arrange: подготавливаем систему и стартовую точку вне ООФ.
    residual_domain_violation_var equation;
    fixed_solver_parameters_t<-1, 0, golden_section_search> parameters;
    Eigen::VectorXd initial(1);
    initial(0) = -1.0;

    fixed_solver_result_t<-1> result;

    // Act: запускаем Ньютон-Рафсон.
    EXPECT_NO_THROW(fixed_newton_raphson<-1>::solve(equation, initial, parameters, &result));

    // Assert: получаем код численной ошибки без проброса исключения наружу.
    ASSERT_EQ(result.result_code, numerical_result_code_t::NumericalNanValues);
}

/// @brief Проверяет перехват domain_violation при расчете Якоби для fixed_system_t<-1>.
TEST(NewtonRaphsonDomainViolation, CatchesJacobianViolationVar)
{
    // Arrange: подготавливаем систему с доменным исключением в jacobian_sparse.
    jacobian_domain_violation_var equation;
    fixed_solver_parameters_t<-1, 0, golden_section_search> parameters;
    Eigen::VectorXd initial(1);
    initial(0) = 0.5;

    fixed_solver_result_t<-1> result;

    // Act: запускаем Ньютон-Рафсон.
    EXPECT_NO_THROW(fixed_newton_raphson<-1>::solve(equation, initial, parameters, &result));

    // Assert: проверяем, что исключение перехвачено внутри алгоритма.
    ASSERT_EQ(result.result_code, numerical_result_code_t::NumericalNanValues);
}

/// @brief Проверяет перехват domain_violation при недопустимом начальном аргументе для fixed_system_t<2>.
TEST(NewtonRaphsonDomainViolation, CatchesResidualViolationAtInitialPointFixed2)
{
    // Arrange: подготавливаем систему и стартовую точку вне ООФ.
    residual_domain_violation_fixed2 equation;
    fixed_solver_parameters_t<2, 0, golden_section_search> parameters;
    std::array<double, 2> initial{ -1.0, 0.0 };

    fixed_solver_result_t<2> result;

    // Act: запускаем Ньютон-Рафсон.
    EXPECT_NO_THROW(fixed_newton_raphson<2>::solve(equation, initial, parameters, &result));

    // Assert: фиксируем текущий контракт по коду результата.
    ASSERT_EQ(result.result_code, numerical_result_code_t::NumericalNanValues);
}

/// @brief Проверяет перехват domain_violation при расчете Якоби для fixed_system_t<2>.
TEST(NewtonRaphsonDomainViolation, CatchesJacobianViolationFixed2)
{
    // Arrange: подготавливаем систему с доменным исключением в jacobian_dense.
    jacobian_domain_violation_fixed2 equation;
    fixed_solver_parameters_t<2, 0, golden_section_search> parameters;
    std::array<double, 2> initial{ 0.0, 0.0 };

    fixed_solver_result_t<2> result;

    // Act: запускаем Ньютон-Рафсон.
    EXPECT_NO_THROW(fixed_newton_raphson<2>::solve(equation, initial, parameters, &result));

    // Assert: проверяем ожидаемый код результата.
    ASSERT_EQ(result.result_code, numerical_result_code_t::NumericalNanValues);
}

/// @brief При исследовании целевой функции в режиме record_nan выход за ООФ отмечается NaN и сохранением индекса
TEST(NewtonRaphsonDomainViolation, RecordsNanAndIndices_WhenLineSearchExploreHitsDomain)
{
    // Arrange: система с границей x<=1 и невязкой с корнем вне ООФ
    residual_domain_violation_var equation;
    equation.upper_bound_for_domain = 1.0;
    equation.linear_residual_root = 1.2;
    fixed_solver_parameters_t<-1, 0, golden_section_search> parameters;
    parameters.iteration_count = 3;
    parameters.analysis.line_search_explore = true;
    parameters.analysis.line_search_explore_on_domain_violation =
        line_search_explore_domain_violation_action_t::record_nan;
    Eigen::VectorXd initial(1);
    initial(0) = 0.95;
    fixed_solver_result_t<-1> result;
    fixed_solver_result_analysis_t<-1> analysis;

    // Act: запускаем Ньютона с аналитикой.
    EXPECT_NO_THROW(fixed_newton_raphson<-1>::solve(equation, initial, parameters, &result, &analysis));

    // Assert: есть согласованные ряды explore и хотя бы один domain_violation на сетке.
    ASSERT_EQ(analysis.target_function.size(), analysis.line_search_explore_domain_violation_indices.size());
    bool found_violation = false;
    for (size_t step = 0; step < analysis.target_function.size(); ++step) {
        const auto& indices = analysis.line_search_explore_domain_violation_indices[step];
        if (indices.empty()) {
            continue;
        }
        found_violation = true;
        for (size_t grid_index : indices) {
            ASSERT_LT(grid_index, analysis.target_function[step].size());
            ASSERT_FALSE(std::isfinite(analysis.target_function[step][grid_index]));
        }
    }
    ASSERT_TRUE(found_violation);
}

/// @brief При исследовании целевой функции в режиме rethrow пробрасывается исключениеdomain_violation
TEST(NewtonRaphsonDomainViolation, PropagatesDomainViolation_WhenLineSearchExploreRethrow)
{
    // Arrange: та же геометрия, режим rethrow.
    residual_domain_violation_var equation;
    equation.upper_bound_for_domain = 1.0;
    equation.linear_residual_root = 1.2;
    fixed_solver_parameters_t<-1, 0, golden_section_search> parameters;
    parameters.iteration_count = 3;
    parameters.analysis.line_search_explore = true;
    parameters.analysis.line_search_explore_on_domain_violation =
        line_search_explore_domain_violation_action_t::rethrow;
    Eigen::VectorXd initial(1);
    initial(0) = 0.95;
    fixed_solver_result_t<-1> result;

    // Act: запускаем Ньютона с исследованием.
    // Assert: ожидаем проброс исключения domain_violation.
    EXPECT_THROW(
        fixed_newton_raphson<-1>::solve(equation, initial, parameters, &result),
        domain_violation
    );
}
