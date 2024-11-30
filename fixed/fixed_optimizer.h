#pragma once


/// @brief Целевая функция суммы квадратов - базовый класс
class fixed_least_squares_function_t {
protected:
    /// @brief Относительное (!) приращение для расчета производных
    double epsilon{ 1e-6 };
public:
    typedef VectorXd argument_type;
    typedef VectorXd residuals_type;
    typedef MatrixXd matrix_type;

public:
    virtual double objective_function(const residuals_type& r) {
        double sum_of_squares = r.squaredNorm();
        return sum_of_squares;
    }
    /// @brief Расчет целевой функции по аргументу
    double operator()(const argument_type& x) {
        residuals_type r = residuals(x);
        return objective_function(r);
    }
    /// @brief Невязки системы уравнений
    virtual residuals_type residuals(const argument_type& x) = 0;
    /// @brief Якобиан системы уравнений
    virtual matrix_type jacobian_dense(const argument_type& x) {
        return jacobian_dense_numeric(x);
    }
    /// @brief Численный расчет якобиана
    matrix_type jacobian_dense_numeric(const argument_type& x)
    {
        argument_type arg = x;

        matrix_type J;

        for (int arg_index = 0; arg_index < x.size(); ++arg_index) {
            double e = epsilon * std::max(1.0, abs(arg[arg_index]));
            arg[arg_index] = x[arg_index] + e;
            residuals_type f_plus = residuals(arg);
            arg[arg_index] = x[arg_index] - e;
            residuals_type f_minus = residuals(arg);
            arg[arg_index] = x[arg_index];

            residuals_type Jcol = (f_plus - f_minus) / (2 * e);

            if (arg_index == 0) {
                J = matrix_type(Jcol.size(), x.size());
            }
            for (size_t row = 0; row < static_cast<size_t>(Jcol.size()); ++row) {
                J(row, arg_index) = Jcol[row];
            }
        }
        return J;
    }
    /// @brief Специфический критерий успешного завершения расчета
    /// @param r Текущее значения невязок
    /// @param x Текущее значение аргумента
    /// @return Флаг успешного завершения
    virtual bool custom_success_criteria(const residuals_type& r, const argument_type& x)
    {
        return true;
    }
};

/// @brief Алгоритм оптимизации Гаусса-Ньютона
/// 2006 Nocedal, Wright 
class fixed_optimize_gauss_newton {
public:
    typedef VectorXd argument_type;
    typedef VectorXd residuals_type;
    typedef MatrixXd matrix_type;
private:
    /// @brief Проверка значения на Nan/infinite для скалярного случая
    /// @param value Проверяемое значение
    /// @return true/false
    static inline bool has_not_finite(const double value)
    {
        if (!std::isfinite(value)) {
            return true;
        }
        return false;
    }

    /// @brief Проверка значения на Nan/infinite для векторного случая
    /// @param value Проверяемое значение
    /// @return true/false
    template <typename Container>
    static inline bool has_not_finite(const Container& values)
    {
        for (const double value : values) {
            if (has_not_finite(value))
                return true;
        }
        return false;
    }
private:
    /// @brief Норма вектора приращения 
    static double argument_increment_factor(
        const argument_type& argument, const argument_type& argument_increment)
    {
        return argument_increment.norm() / argument.size();

    }
private:
    /// @brief Проведение процедуры линейного поиска по заданному алгоритму
    /// @tparam LineSearch Алгоритм линейного поиска
    /// @param line_search_parameters параметры линейного поиска
    /// @param residuals Невязки
    /// @param argument Текущий аргумент, относительного которого делается приращение
    /// @param r Текущая невязка в системе уравнений для быстрого расчета ц.ф.
    /// @param p Приращение аргумента
    /// @return Величина шага. По соглашению если алгоритм линейного поиска не сошелся, будет NaN
    template <typename LineSearch>
    static double perform_line_search(
        const typename LineSearch::parameters_type& line_search_parameters,
        fixed_least_squares_function_t& function,
        const argument_type& argument, const residuals_type& r, const argument_type& p)
    {
        auto directed_function = [&](double step) {
            return function(argument + step * p);
            };

        // Диапазон поиска, значения функции на границах диапазона
        double a = 0;
        double b = line_search_parameters.maximum_step;
        double function_a = directed_function(a);
        double function_b = directed_function(b);

        auto [search_step, elapsed_iterations] = LineSearch::search(
            line_search_parameters,
            directed_function, a, b, function_a, function_b);
        return search_step;
    }

public:
    /// @brief Запуск численного метода
    /// @tparam LineSearch Алгоритм регулировки шага поиска
    /// @param residuals Функция невязок
    /// @param initial_argument Начальное приближение
    /// @param solver_parameters Настройки поиска
    /// @param result Результаты расчета
    template <
        typename LineSearch = divider_search>
    static void optimize(
        fixed_least_squares_function_t& function,
        const argument_type& initial_argument,
        const fixed_solver_parameters_t<-1, 0, LineSearch>& solver_parameters,
        fixed_solver_result_t<-1>* result,
        fixed_solver_result_analysis_t<-1>* analysis = nullptr
    )
    {
        // Информация о том, как изменялась целевая функция по траектории шага
        if (analysis != nullptr && solver_parameters.analysis.line_search_explore) {
            throw std::runtime_error("Line search exploring not implemented");
            //analysis->target_function.push_back(perform_step_research(residuals, argument, p));
        }

        size_t n = initial_argument.size();

        if (n == 0) {
            result->result_code = numerical_result_code_t::Converged;
            return;
        }

        result->argument = initial_argument;
        VectorXd& argument = result->argument;
        VectorXd& r = result->residuals;

        if (analysis != nullptr && solver_parameters.analysis.argument_history) {
            analysis->argument_history.push_back(argument);
        }

        VectorXd argument_increment, search_direction;

        r = function.residuals(argument);
        result->residuals_norm = function.objective_function(r);
        if (analysis != nullptr && solver_parameters.analysis.objective_function_history) {
            analysis->target_function.push_back({ result->residuals_norm });
        }

        bool custom_stop_criteria = function.custom_success_criteria(r, argument);
        /*if (custom_stop_criteria) {
            result->result_code = numerical_result_code_t::Converged;
            return;
        }*/

        for (size_t iteration = 0; iteration < solver_parameters.iteration_count; ++iteration)
        {
            MatrixXd J = function.jacobian_dense(argument);

            /*       if (solver_parameters.optimize_quadprog) {
                       // для сравнения
                       JacobiSVD<MatrixXd> svd_solver(J, ComputeThinU | ComputeThinV);
                       VectorXd search_direction_unconstrained = svd_solver.solve(-r);
                       callback.after_search_direction(iteration, &argument, &r, &Jsparse, &search_direction_unconstrained);

                       MatrixXd H = J.transpose() * J;
                       VectorXd f = r.transpose() * J;

                       MatrixXd A;  VectorXd b;
                       solver_parameters.constraints.get_inequalities_constraints(argument, &A, &b);

                       // см. формулы в "Идентификация - quadprog.docx"
                       MatrixXd A_increment = A;
                       VectorXd b_increment = -A * argument + b;
                       solve_quadprog_inequalities(H, f, A_increment, b_increment,
                           VectorXd::Zero(argument.size()), &search_direction);

                       trim_increment(solver_parameters.constraints.boundaries, search_direction);
                       trim_increment_relative(solver_parameters.constraints.relative_boundaries, search_direction);


                       callback.after_search_direction_trim(iteration, &argument, &r, &Jsparse, &search_direction);
                   }
                   else */

            {
                // Наш Якобиан J получен из ряда Тейлора r(x0 + dx) = r(x0) + J(x0)dx, 
                // линеаризованная задача МНК выглядит так: ||r(x0) + J(x0)dx|| -> min по dx
                // Алгоритм JacobiSVD подразумевает задачу  ||Y - Ax|| -> min по x
                // В итоге, либо передавать -J, либо -r. Последнее вычислительно быстрее
                Eigen::JacobiSVD<MatrixXd> svd_solver(J, Eigen::ComputeThinU | Eigen::ComputeThinV);
                search_direction = svd_solver.solve(-r);

                //trim_increment_min(solver_parameters.constraints.minimum, argument, search_direction);
                //trim_increment_max(solver_parameters.constraints.maximum, argument, search_direction);
                //trim_increment(solver_parameters.constraints.boundaries, search_direction);
                //trim_increment_relative(solver_parameters.constraints.relative_boundaries, search_direction);


            }


            double search_step = perform_line_search<LineSearch>(
                solver_parameters.line_search, function, argument, r, search_direction);
            if (analysis != nullptr && solver_parameters.analysis.steps) {
                analysis->steps.push_back(search_step);
            }

            if (!std::isfinite(search_step)) {
                result->result_code = numerical_result_code_t::LineSearchFailed;
                break;
            }

            argument_increment = search_step * search_direction;
            argument += argument_increment;

            if (analysis != nullptr && solver_parameters.analysis.argument_history) {
                analysis->argument_history.push_back(argument);
            }

            r = function.residuals(argument);
            result->residuals_norm = function.objective_function(r);
            if (analysis != nullptr && solver_parameters.analysis.objective_function_history) {
                analysis->target_function.push_back({ result->residuals_norm });
            }

            // Проверка критерия выхода по малому относительному приращению
            double argument_increment_metric = solver_parameters.step_criteria_assuming_search_step
                ? argument_increment_factor(argument, argument_increment)
                : argument_increment_factor(argument, search_direction);
            bool argument_increment_criteria =
                argument_increment_metric < solver_parameters.argument_increment_norm;
            bool custom_criteria = false;
            if (custom_criteria || argument_increment_criteria) {
                result->result_code = numerical_result_code_t::Converged;
                break;
            }

            //double argument_increment_metric = solver_parameters.step_criteria_assuming_search_step
            //    ? argument_increment_factor(gains, argument_increment)
            //    : argument_increment_factor(gains, search_direction);

            //result->argument_increment_criteria =
            //    argument_increment_metric < solver_parameters.argument_increment_norm;
            //result->custom_stop_criteria = function.custom_stop_criteria(r);
            //result->enough_iterations_criteria = iteration + 1 == solver_parameters.iteration_count;

            //if (result->has_any_stop_criteria()) {
            //    result->converged = true;
            //    break;
            //}
        }
    };

};

/// @brief Упрощенный запуск оптимизатора Гаусса-Ньютона
/// 
/// @param function Целевая функция
/// @param initial Начальное приближение
/// @return Результат оптимизации
template <typename ResidualFunction>
inline VectorXd optimize_gauss_newton(ResidualFunction& function, const VectorXd& initial)
{
    fixed_solver_parameters_t<-1, 0, golden_section_search> parameters;
    fixed_solver_result_t<-1> result;

    fixed_optimize_gauss_newton::optimize(function, initial, parameters, &result);
    if (result.result_code != numerical_result_code_t::Converged) {
        throw std::runtime_error("Gauss-Newton optimizer not converged");
    }

    return result.argument;
};
