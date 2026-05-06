#pragma once

namespace {
/// @brief Система переменной размерности, кидающая domain_violation при недопустимом x.
struct residual_domain_violation_var : public fixed_system_t<-1> {
    Eigen::VectorXd residuals(const Eigen::VectorXd& x) override {
        if (x(0) < 0.0) {
            throw domain_violation{};
        }
        Eigen::VectorXd r(1);
        r(0) = x(0) - 1.0;
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
