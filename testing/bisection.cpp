#include "gtest/gtest.h"
#include <fixed/fixed_bisection.h>


/// \brief линейное уравнение
struct simple_linear : public fixed_system_t<1> {
public:
    /// \brief
    /// функция вычисления невязки
    /// \param x - аргумент
    /// \return невязка
    virtual function_type residuals(const var_type &x) override
    {
        return 5.-x;
    }
    /// \brief ожидаемое решение
    /// \return ожидаемое решение
    constexpr fixed_system_t<1>::var_type expectet_solution() const { return 5.;}
};

/// \brief квадратное уравнение - пример для теста
struct simple_quadratic : public fixed_system_t<1> {
public:
    /// \brief расчет невязки квадратного уравнения
    virtual function_type residuals(const var_type &x) override
    {
        return 3.5-x*x;
    }
    /// \brief ожидаемое решение
    /// \return ожидаемое решение
    const fixed_system_t<1>::var_type expectet_solution() const { return std::sqrt(3.5);}
};

/// \brief
/// тест проверки работоспособности метода бисекции в случае
/// линейного уравнения
///
TEST(Bisection, SolvesLinearEquation){
    simple_linear eqn;
    fixed_bisectional_parameters_t p;
    p.argument_limit_min=0;
    p.argument_limit_max=7.8;
    p.argument_history = true;
    p.secant_treshhold_iterations = 0;
    p.solution_type = fixed_bisectional_solution_type::Secant;
    p.verbose=false;
    fixed_bisection_result_t<1> res;
    fixed_bisection_result_analysis_t<1> ana;
    fixed_bisectional<1>::solve(p, eqn, &res, &ana);

    EXPECT_DOUBLE_EQ(res.argument,eqn.expectet_solution());
    if(p.verbose){
        std::cout<<res.argument<<'\t'<<res.result_code<<'\t'<<res.score<<'\t'<<res.iteration_count<<std::endl;
        for(auto a:ana.argument_history){
            std::cout<<'\t'<<a<<std::endl;
        }
    }
};

/// \brief
/// тест проверки работоспособности метода бисекции в случае
/// квадратного уравнения
TEST(Bisection, SolvesQuadraticEquation){
    simple_quadratic eqn;
    fixed_bisectional_parameters_t p;
    p.argument_limit_min=0;
    p.argument_limit_max=27.8;
    p.argument_history = true;
    p.solution_type = fixed_bisectional_solution_type::Combined;
    p.secant_treshhold_max =0.01;
    p.secant_treshhold_min =0.001;
    p.verbose=false;
    fixed_bisection_result_t<1> res;
    fixed_bisection_result_analysis_t<1> ana;
    fixed_bisectional<1>::solve(p, eqn, &res, &ana);

    EXPECT_NEAR(res.argument,eqn.expectet_solution(),p.residual_precision);
    if(p.verbose){
        std::cout<<res.argument<<'\t'<<res.result_code<<'\t'<<res.score<<'\t'<<res.iteration_count<<std::endl;
        for(auto a:ana.argument_history){
            std::cout<<'\t'<<a<<std::endl;
        }
    }
};

