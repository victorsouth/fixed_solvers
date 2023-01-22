#pragma once

#include "line_search/golden_section.h"
#include "line_search/divider.h"


/// @brief Ограничения солвера
template <std::ptrdiff_t Dimension>
struct fixed_solver_constraints;

/// @brief Ограничения солвера для переменной размерности 
/// (пока не используется, заготовка на будущее)
template <>
struct fixed_solver_constraints<-1>
{
    std::vector<std::pair<size_t, double>> relative_boundary;
    std::vector<std::pair<size_t, double>> minimum;
    std::vector<std::pair<size_t, double>> maximum;
    size_t get_constraint_count() const
    {
        return minimum.size() + maximum.size();
    }

    static std::pair<MatrixXd, VectorXd> get_inequalities_constraints_vectors(
        size_t argument_dimension, 
        const std::vector<std::pair<size_t, double>>& boundaries)
    {
        MatrixXd A = MatrixXd::Zero(boundaries.size(), argument_dimension);
        VectorXd b = VectorXd::Zero(boundaries.size());

        int row_index = 0;
        for (const auto& kvp : boundaries) {
            A(row_index, kvp.first) = 1;
            b(row_index) = kvp.second;
            row_index++;
        }
        return std::make_pair(A, b);
    }

    /// @brief Учитывает только minimum, maximum
    std::pair<MatrixXd, VectorXd> get_inequalities_constraints(const size_t argument_size) const
    {
        // ограничения
        const auto n = argument_size;
        MatrixXd A = MatrixXd::Zero(get_constraint_count(), n);
        VectorXd B = VectorXd::Zero(get_constraint_count());

        size_t offset = 0;
        // максимальное значение
        { 
            auto [a, b] = get_inequalities_constraints_vectors(n, maximum);

            A.block(offset, 0, a.rows(), n) = a;
            B.segment(offset, b.size()) = b;
            offset += a.rows();
        }
        // минимальное значение
        {
            auto [a, b] = get_inequalities_constraints_vectors(n, minimum);

            A.block(offset, 0, a.rows(), n) = -a;
            B.segment(offset, b.size()) = -b;
            offset += a.rows();
        }
        return std::make_pair(A, B);
    }


    void trim_increment_relative(VectorXd& increment) const
    {
        double factor = 1;
        for (const auto& kvp : relative_boundary)
        {
            int sign = sgn(increment(kvp.first));
            if (sign * increment(kvp.first) > kvp.second) {
                double current_factor = sign * increment(kvp.first) / kvp.second;
                factor = std::max(factor, current_factor);
            }
        }
        if (factor > 1)
        {
            increment /= factor;
        }
    }

};


/// @brief Ограничения солвера
template <std::ptrdiff_t Dimension>
struct fixed_solver_constraints
{
    typedef typename fixed_system_types<Dimension>::var_type var_type;
    var_type relative_boundary{ fixed_system_types<Dimension>::default_var() };
    var_type minimum{ fixed_system_types<Dimension>::default_var() };
    var_type maximum{ fixed_system_types<Dimension>::default_var() };

    /// @brief Обрезка по минимуму для скалярного случая
    void trim_min(double argument, double& increment) const
    {
        if (std::isnan(minimum))
            return;

        if (argument + increment < minimum) {
            increment = minimum - argument;
        }
    }
    /// @brief Обрезка по минимуму для векторного случая
    void trim_min(const array<double, Dimension>& argument,
        array<double, Dimension>& increment) const
    {
        double factor = 1;

        for (size_t index = 0; index < increment.size(); ++index) {
            if (std::isnan(minimum[index])) {
                continue;
            }

            if (argument[index] + increment[index] < minimum[index]) {
                double allowed_increment = minimum[index] - argument[index];
                factor = std::max(factor, increment[index] / allowed_increment);
            }
        }
        if (factor > 1)
        {
            increment = increment / factor;
        }
    }

    /// @brief Приводит значение аргумента внутрь ограничений мин/макс
    void ensure_constraints(double& argument) const
    {
        if (!std::isnan(maximum)) {
            argument = std::min(argument, maximum);
        }
        if (!std::isnan(minimum)) {
            argument = std::max(argument, minimum);
        }
    }


    /// @brief Обрезка по минимуму для скалярного случая
    void trim_max(double argument, double& increment) const
    {
        if (std::isnan(maximum))
            return;

        if (argument + increment > maximum) {
            increment = maximum - argument;
        }
    }

    void trim_max(const array<double, Dimension>& argument,
        array<double, Dimension>& increment) const
    {
        double factor = 1;

        for (size_t index = 0; index < increment.size(); ++index) {
            if (std::isnan(minimum[index])) {
                continue;
            }

            if (argument[index] + increment[index] > maximum[index]) {
                double allowed_increment = maximum[index] - argument[index];
                factor = std::max(factor, increment[index] / allowed_increment);
            }
        }
        if (factor > 1)
        {
            increment = increment / factor;
        }
    }

    /// @brief Обрезка по максимального приращению для скларяного случая
    void trim_relative(double& increment) const
    {
        if (std::isnan(relative_boundary))
            return;
        double abs_inc = abs(increment);
        if (abs_inc > relative_boundary) {
            double factor = abs_inc / relative_boundary;
            increment /= factor;
        }
    }
    /// @brief Обрезка по максимального приращению для векторного случая
    void trim_relative(array<double, Dimension>& increment) const
    {
        double factor = 1;
        for (size_t index = 0; index < increment.size(); ++index) {
            if (std::isnan(relative_boundary[index])) {
                continue;
            }

            double abs_inc = abs(increment[index]);
            if (abs_inc > relative_boundary[index]) {
                double current_factor = abs_inc / relative_boundary[index];
                factor = std::max(factor, current_factor);
            }
        }
        if (factor > 1)
        {
            increment = increment / factor;
        }
    }
};

/// @brief Параметры диагностики сходимости
struct fixed_solver_analysis_parameters_t {
    /// @brief Собирает историю значений аргумента
    bool argument_history{ false };
    /// @brief Собирать значения шага в методе Ньютона-Рафсона
    bool steps{false};
    /// @brief Выполнять исследование целевой функции на каждом шаге
    bool line_search_explore{ false };
};

/// @brief Параметры алгоритма Ньютона-Рафсона
template <std::ptrdiff_t Dimension, typename LineSearch = divider_search>
struct fixed_solver_parameters_t
{
    /// @brief Параметры диагностики
    fixed_solver_analysis_parameters_t analysis;
    /// @brief Границы на диапазон и единичный шаг
    fixed_solver_constraints<Dimension> constraints;
    /// @brief Параметры алгоритма регулировки шага
    typename LineSearch::parameters_type line_search;
    /// @brief Количество итераций
    size_t iteration_count{ 100 };
    /// @brief Погрешность метода по приращению аргумента
    double argument_increment_norm{ 1e-4 };
    /// @brief Проверять векторный шага после регулировки или перед ней
    bool step_criteria_assuming_search_step{ false };
    /// @brief Считать невозможность улучшить ц.ф. при линейном поиске как лучшее решений
    /// Т.е. считать, что причина этого в наличии шума уравнений и ц.ф.
    bool treat_line_search_fail_as_convergence{ false };
};

enum class numerical_result_code_t
{
    NoNumericalError, IllConditionedMatrix, LargeConditionNumber,
    NotConverged, NumericalNanValues, LineSearchFailed, Converged
};

/// @brief Результат расчета численного метода
template <std::ptrdiff_t Dimension>
struct fixed_solver_result_t {
    typedef typename fixed_system_types<Dimension>::var_type var_type;
    typedef typename fixed_system_types<Dimension>::right_party_type function_type;

    /// @brief Метрика приращения
    double argument_increment_metric{ 0 };
    /// @brief Показывает, что достигнуто минимальное приращение спуска
    bool argument_increment_criteria{ false };
    /// @brief Завершение итерационной процедуры
    numerical_result_code_t result_code{ numerical_result_code_t::NotConverged };
    /// @brief Остаточная невязка по окончании численного метода
    function_type residuals;
    /// @brief Искомый аргумент по окончании численного метода
    var_type argument;
    /// @brief Количество итераций
    size_t iteration_count{ 0 };
};


/// @brief Аналитика сходимости
template <size_t Dimension>
struct fixed_solver_result_analysis_t {
public:
    /// @brief Значения целевой функции для одной регулировки шага
    typedef vector<double> target_function_values_t;
    typedef typename fixed_system_types<Dimension>::var_type var_type;
    typedef typename fixed_system_types<Dimension>::right_party_type function_type;
    typedef typename fixed_system_types<Dimension>::equation_coeffs_type equation_coeffs_type;
public:
    /// @brief Значения целевой функции по всем шагам Ньютона-Рафсона
    vector<target_function_values_t> target_function;
    /// @brief Результат расчета
    vector<var_type> argument_history;
    /// @brief Величина шагов Ньютона-Рафсона
    vector<double> steps;
};


/// @brief Солвер Ньютона-Рафсона для систем уравнений фиксированной размерности
/// В том числе скалярных
/// Ньютона-Рафсона - означает регулировку шага
/// Реализована регулировка шага за счет дробления
template <size_t Dimension>
class fixed_newton_raphson {
public:
    typedef typename fixed_system_types<Dimension>::var_type var_type;
    typedef typename fixed_system_types<Dimension>::right_party_type function_type;
    typedef typename fixed_system_types<Dimension>::equation_coeffs_type equation_coeffs_type;
private:
    static inline bool has_not_finite(const double value)
    {
        if (!std::isfinite(value)) {
            return true;
        }
        return false;
    }

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
    static double argument_increment_factor(
        const var_type& argument, const var_type& argument_increment);
private:
    /// @brief Расчет целевой функции при регулировке шага в диапазоне [0, 1]
    /// @param residuals Векторная функция невязок 
    /// @param argument Текущее значение аргумента
    /// @param p Значение приращения по методу Ньютона
    /// @return Результат исследования ц.ф. 
    static vector<double> perform_step_research(
        fixed_system_t<Dimension>& residuals,
        const var_type& argument,
        const var_type& p)
    {
        auto directed_function = [&](double step) {
            return residuals(argument + step * p);
        };

        size_t research_step_count = 100;
        vector<double> target_function;
        target_function.reserve(research_step_count);
        for (size_t index = 0; index <= research_step_count; ++index) {
            double alpha = 1.0 * index / research_step_count;
            double norm = directed_function(alpha);
            target_function.emplace_back(norm);
        }

        return std::move(target_function);
    }


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
        fixed_system_t<Dimension>& residuals, 
        const var_type& argument, const var_type& r, const var_type& p) 
    {
        auto directed_function = [&](double step) {
            return residuals(argument + step * p);
        };

        // Диапазон поиска, значения функции на границах диапазона
        double a = 0;
        double b = line_search_parameters.maximum_step;
        double function_a = residuals.objective_function(r); // уже есть рассчитанные невязки, по ним считаем ц.ф.
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
    template <typename LineSearch = divider_search>
    static void solve_dense(
        fixed_system_t<Dimension>& residuals,
        const var_type& initial_argument,
        const fixed_solver_parameters_t<Dimension, LineSearch>& solver_parameters,
        fixed_solver_result_t<Dimension>* result, 
        fixed_solver_result_analysis_t<Dimension>* analysis = nullptr
    )
    {
        var_type& r = result->residuals;
        function_type& argument = result->argument;
        double& argument_increment_metric = result->argument_increment_metric;
        bool& argument_increment_criteria = result->argument_increment_criteria;


        argument = initial_argument;
        if (analysis != nullptr && solver_parameters.analysis.argument_history) {
            analysis->argument_history.push_back(argument);
        }

        var_type argument_increment;
        var_type p;

        r = residuals.residuals(argument);
        if (has_not_finite(r)) {
            r = residuals.residuals(argument); // для отладки
            result->result_code = numerical_result_code_t::NumericalNanValues;
            return;
        }

        result->result_code = numerical_result_code_t::NotConverged;

        size_t& iteration = result->iteration_count;
        for (iteration = 0; iteration < solver_parameters.iteration_count; ++iteration)
        {
            auto J = residuals.jacobian_dense(argument);
            p = -solve_linear_system(J, r);
            // todo: обработчик ошибки решения СЛАУ

            if (solver_parameters.step_criteria_assuming_search_step == false)
            {
                argument_increment_metric = argument_increment_factor(argument, p);
                argument_increment_criteria =
                    argument_increment_metric < solver_parameters.argument_increment_norm;
                if (argument_increment_criteria) {
                    result->result_code = numerical_result_code_t::Converged;
                    break;
                }
            }

            solver_parameters.constraints.trim_max(argument, p);
            solver_parameters.constraints.trim_min(argument, p);
            if constexpr (std::is_same<LineSearch, divider_search>()) {
                // Для ЗС trim делать не надо
                solver_parameters.constraints.trim_relative(p);
            }

            if (analysis != nullptr && solver_parameters.analysis.line_search_explore) {
                analysis->target_function.push_back(perform_step_research(residuals, argument, p));
            }

            double search_step = perform_line_search<LineSearch>(
                solver_parameters.line_search, residuals, argument, r, p);
            if (std::isnan(search_step)) {
                if (solver_parameters.treat_line_search_fail_as_convergence) {
                    result->result_code = numerical_result_code_t::Converged;
                }
                else {
                    result->result_code = numerical_result_code_t::LineSearchFailed;
                }
                break;
            }

            if (analysis != nullptr && solver_parameters.analysis.steps) {
                analysis->steps.push_back(search_step);
            }
            

            argument_increment = search_step * p;
            argument += argument_increment;

            if (analysis != nullptr && solver_parameters.analysis.argument_history) {
                analysis->argument_history.push_back(argument);
            }

            r = residuals.residuals(argument);
            if (has_not_finite(r)) {
                r = residuals.residuals(argument); // для отладки
                result->result_code = numerical_result_code_t::NumericalNanValues;
                break;
            }

            argument_increment_metric = solver_parameters.step_criteria_assuming_search_step
                ? argument_increment_factor(argument, argument_increment)
                : argument_increment_factor(argument, p);
            argument_increment_criteria =
                argument_increment_metric < solver_parameters.argument_increment_norm;
            bool custom_criteria = false;
            if (custom_criteria || argument_increment_criteria) {
                result->result_code = numerical_result_code_t::Converged;
                break;
            }
        }
    }

};

template <>
inline double fixed_newton_raphson<1>::argument_increment_factor(
    const double& argument,
    const double& argument_increment)
{
    double arg = std::max(1.0, argument);
    double inc = argument_increment;

    double result = abs(inc / arg);
    return result;
}

template <size_t Dimension>
inline double fixed_newton_raphson<Dimension>::argument_increment_factor(
    const var_type& argument, const var_type& argument_increment)
{
    double squared_sum = 0;

    for (int component = 0; component < argument_increment.size(); ++component) {
        double arg = std::max(1.0, argument[component]);
        double inc = argument_increment[component];
        squared_sum += pow(inc / arg, 2);
    }
    return sqrt(squared_sum) / Dimension;
}
