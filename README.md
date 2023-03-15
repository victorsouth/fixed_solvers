# fixed_solvers
Данный солвер решает систему нелинейных уравнений фикисрованной размерности с помощью алгоритма Ньютона-Рафсона.
Учет фиксированности размерности позволяет не использовать динамическую память для хранения векторов.
---------------------------------------------------------------------------------------------------------------
class sample_system : public fixed_system_t<2>
    {
        using fixed_system_t<2>::var_type;
    public:
        var_type residuals(const var_type& x) {
            return
            {
                pow(x[0] - 2, 3),
                pow(x[1] - 1, 3)
            };
        }
    };

    sample_system test;
    fixed_solver_parameters_t<2, 0> parameters;
    fixed_solver_result_t<2> result;
    fixed_newton_raphson<2>::solve_dense(test, { 0, 0 }, parameters, & result);
---------------------------------------------------------------------------------------------------------------
