﻿#pragma once

#include "line_search/golden_section.h"
#include "line_search/divider.h"

using std::vector;
using std::map;
using std::wstringstream;
using std::wstring;

/// @brief Точность проверки нахождения на ограничениях
constexpr double eps_constraints = 1e-8;

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

                
                if (std::abs(argument[index] - minimum[index]) < eps_constraints) {
                    // Параметр уже сел на ограничения, allowed_decrement будет нулевой,
                    // соответственно factor получается бесконечный (см. ветвь "else" ниже)
                    // не учитываем эту переменную при расчете factor, сразу обрезаем 
                    increment[index] = 0;
                }
                else {
                    double allowed_decrement = argument[index] - minimum[index];
                    factor = std::max(factor, abs(increment[index]) / allowed_decrement);
                }

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
            if (std::isnan(maximum[index])) {
                continue;
            }

            if (argument[index] + increment[index] > maximum[index]) {
                if (std::abs(argument[index] - maximum[index]) < eps_constraints) {
                    // Параметр уже сел на ограничения, allowed_increment будет нулевой,
                    // соответственно factor получается бесконечный (см. ветвь "else" ниже)
                    // не учитываем эту переменную при расчете factor, сразу обрезаем 
                    increment[index] = 0;
                }
                else {
                    double allowed_increment = maximum[index] - argument[index];
                    factor = std::max(factor, increment[index] / allowed_increment);
                }
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

/// @brief Действия при неудачной регулировке шага
enum class line_search_fail_action_t { TreatAsFail, TreatAsSuccess, PerformMinStep };

/// @brief Линейные ограничения - произвольная размерность
template <std::ptrdiff_t Dimension, std::ptrdiff_t Count>
struct fixed_linear_constraints;

/// @brief Специализация линейных ограничений 
/// для случая отсутствия линейныйх ограничений (их ноль)
template <std::ptrdiff_t Dimension>
struct fixed_linear_constraints<Dimension, 0> {
    /// @brief Функция trim ничего не делает, ее просто можно вызвать
    typedef typename fixed_system_types<Dimension>::var_type var_type;
    inline void trim(const var_type& argument, var_type& increment) const
    { }
};

/// @brief Специализация для систем второго порядка, одно ограничение
template <>
struct fixed_linear_constraints<2, 1> {
    typedef typename fixed_system_types<2>::var_type var_type;
    typedef typename fixed_system_types<2>::matrix_type matrix_type;
    /// @brief Коэффициенты a (левая часть ax <= b)
    var_type a{ std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN() };
    /// @brief Коэффициенты b (правой часть ax <= b)
    double b{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief Проверяет, что приближение x не нарушает ограничения
    bool check_constraint_satisfaction(const var_type& x) const {
        if (std::isfinite(b))
            return inner_prod(a, x) <= b;
        else
            return true;
    }
    /// @brief Проверяет, что приближение x находится на границе ограничений
    bool check_constraint_border(const var_type& x) const {
        if (std::isfinite(b))
            return std::abs(inner_prod(a, x) - b) < eps_constraints;
        else
            return true;
    }

    /// @brief На основе канонического уравления прямой отдает коэффициенты уравнения
    /// ax = b
    /// @param p1 Первая точка
    /// @param p2  Вторая точка
    /// @return Пара вектор, скаляр: (a, b)
    static std::pair<var_type, double> get_line_coeffs(const var_type& p1, const var_type& p2) {
        double x1 = p1[0];
        double y1 = p1[1];
        double x2 = p2[0];
        double y2 = p2[1];

        // Это в виде y = kx + b
        double k = (y2 - y1) / (x2 - x1);
        double b = y1 - k * x1;

        // -k*x + 1y = b

        return make_pair(var_type{ -k, 1.0 }, b);
    }

    /// @brief Режет по линейным ограниченями a'x <= b
    /// @param argument Текущее значение
    /// @param increment Приращение
    void trim(const var_type& x, var_type& dx) const
    {
        if (!std::isfinite(b))
            return;

        const var_type& p1 = x;
        const var_type p2 = x + dx;

        if (check_constraint_satisfaction(p2))
            return;

        if (check_constraint_border(p1)) {
            // Тут что-то понять можно только по рисунку
            double k = -a[0] / a[1];
            double alpha = atan(k);
            double p = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);
            double beta = acos(dx[0] / p);

            double gamma = beta - alpha;
            double p_dash = p * cos(gamma);
            double px = p_dash * cos(alpha);
            double py = p_dash * sin(alpha);
            dx = { px, py };
        }
        else {
            auto [a2, b2] = get_line_coeffs(p1, p2);
            var_type x_star = solve_linear_system({ a, a2 }, { b, b2 });
            dx = x_star - x;
        }
    }
};

/// @brief Заглушка для ограничений, не даст собраться методу Ньютона, 
/// т.к. отсутствует метод trim
template <std::ptrdiff_t Dimension, std::ptrdiff_t Count>
struct fixed_linear_constraints
{
};


/// @brief Параметры алгоритма Ньютона-Рафсона
template <
    std::ptrdiff_t Dimension, 
    std::ptrdiff_t LinearConstraintsCount = 0, 
    typename LineSearch = divider_search>
struct fixed_solver_parameters_t
{
    /// @brief Параметры диагностики
    fixed_solver_analysis_parameters_t analysis;
    /// @brief Границы на диапазон и единичный шаг
    fixed_solver_constraints<Dimension> constraints;
    /// @brief Линейные ограничения
    fixed_linear_constraints<Dimension, LinearConstraintsCount> linear_constraints;
    /// @brief Параметры алгоритма регулировки шага
    typename LineSearch::parameters_type line_search;
    /// @brief Количество итераций
    size_t iteration_count{ 100 };
    /// @brief Критерий выхода из процедуры (погрешность метода) по приращению аргумента 
    double argument_increment_norm{ 1e-4 };
    /// @brief Критерий выхода из процедуры (погрешность метода) по норме невязок
    double residuals_norm{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief Проверять векторный шага после регулировки или перед ней
    bool step_criteria_assuming_search_step{ false };
    /// @brief Действие при ошибке регулировки шага
    line_search_fail_action_t line_search_fail_action{ line_search_fail_action_t::TreatAsFail};
};

enum class numerical_result_code_t
{
    NoNumericalError, IllConditionedMatrix, LargeConditionNumber, CustomCriteriaFailed,
    NotConverged, NumericalNanValues, LineSearchFailed, Converged
};

/// @brief Оценка качества сходимости
enum class convergence_score_t : int {
    Excellent = 5, Good = 4, Satisfactory = 3, Poor = 2, Error = 1
};

/// @brief Общее количество расчетов при нагрузочном расчете
/// @param score Статистика расчетов. Ключ - оценка, Значение - количество таких полученных оценок
/// @return Общее количество расчетов
inline size_t get_score_total_calculations(const std::map<convergence_score_t, size_t>& score)
{
    size_t result = 0;
    for (auto [sc, count] : score) {
        result += count;
    }
    return result;
}

/// @brief Строковое представление статистики массового расчета
/// В строку выводятся только те оценки, которые реально встретились в статистике
/// @param score Статистика расчетов. Ключ - оценка, Значение - количество таких полученных оценок
/// @return Строка с оценками
inline wstring get_score_string(
    const map<convergence_score_t, size_t>& score)
{
    static const map<convergence_score_t, wstring> score_strings{
        {convergence_score_t::Excellent, L"Exc"},
        {convergence_score_t::Good , L"Good"},
        {convergence_score_t::Satisfactory, L"Satis"},
        {convergence_score_t::Poor, L"Poor"},
        {convergence_score_t::Error, L"Error"},
    };

    size_t total_calculations = get_score_total_calculations(score);

    wstringstream result;
    result << std::setprecision(4);

    for (auto [sc, count] : score) {
        double percent = 100 * ((double)count) / total_calculations;
        result << L"" << score_strings.at(sc) << L": " << percent << L"% ";
    }

    return result.str();
}

/// @brief Процент расчетов, которые успешно завершились 
/// Суммарное количество оценок Отл, Хор, Удовл
/// @param score Статистика расчетов. Ключ - оценка, Значение - количество таких полученных оценок
/// @return Процент (именно процент, не доля) успешных расчетов 
inline double get_converged_percent(
    const map<convergence_score_t, size_t>& score)
{
    size_t total_calculations = get_score_total_calculations(score);

    size_t converged_count = 0;
    for (auto [sc, count] : score) {
        if (sc == convergence_score_t::Excellent ||
            sc == convergence_score_t::Good ||
            sc == convergence_score_t::Satisfactory
            )
        {
            converged_count += count;
        }
    }

    return (double)converged_count / total_calculations * 100;
}


/// @brief Граница малого шага между "отлично" и "хорошо"
constexpr double small_step_threshold{ 0.1 };


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
    /// @brief Балл сходимости
    convergence_score_t score{convergence_score_t::Poor};
    /// @brief Остаточная невязка по окончании численного метода
    function_type residuals{ fixed_system_types<Dimension>::default_var() };
    /// @brief Норма остаточной невязки по окончании численного метода
    double residuals_norm{0};
    /// @brief Искомый аргумент по окончании численного метода
    var_type argument{ fixed_system_types<Dimension>::default_var()};
    /// @brief Количество выполненных итераций
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
    template <
        std::ptrdiff_t LinearConstraintsCount,
        typename LineSearch = divider_search>
    static void solve_dense(
        fixed_system_t<Dimension>& residuals,
        const var_type& initial_argument,
        const fixed_solver_parameters_t<Dimension, LinearConstraintsCount, LineSearch>& solver_parameters,
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
            result->score = convergence_score_t::Error;
            return;
        }

        result->result_code = numerical_result_code_t::NotConverged;
        result->score = convergence_score_t::Excellent;

        size_t& iteration = result->iteration_count;
        for (iteration = 0; iteration < solver_parameters.iteration_count; ++iteration)
        {
            if (iteration > 0.3 * solver_parameters.iteration_count) {
                result->score = std::min(result->score, convergence_score_t::Satisfactory);
            }
            else if (iteration > 0.15 * solver_parameters.iteration_count)
            {
                result->score = std::min(result->score, convergence_score_t::Good);
            }

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

            solver_parameters.linear_constraints.trim(argument, p);
            solver_parameters.constraints.trim_max(argument, p);
            solver_parameters.constraints.trim_min(argument, p);
            solver_parameters.constraints.trim_relative(p);

            if (analysis != nullptr && solver_parameters.analysis.line_search_explore) {
                analysis->target_function.push_back(perform_step_research(residuals, argument, p));
            }

            double search_step = perform_line_search<LineSearch>(
                solver_parameters.line_search, residuals, argument, r, p);
            if (std::isfinite(search_step))
            {
                if (search_step < small_step_threshold) {
                    // снижаем до 4-х баллов, если он выше
                    result->score = std::min(result->score, convergence_score_t::Good);
                }
            }
            else
            {
                // Обработка невозможности спуска ц.ф. при линейном поиске
                if (solver_parameters.line_search_fail_action == line_search_fail_action_t::PerformMinStep) {
                    // понижаем балл до удовлетворительно, считаем дальше
                    result->score = std::min(result->score, convergence_score_t::Satisfactory);
                    search_step = solver_parameters.line_search.step_on_search_fail();
                }
                else if (solver_parameters.line_search_fail_action == line_search_fail_action_t::TreatAsFail) {
                    result->result_code = numerical_result_code_t::LineSearchFailed;
                    break;
                }
                else if (solver_parameters.line_search_fail_action == line_search_fail_action_t::TreatAsSuccess) {
                    result->result_code = numerical_result_code_t::Converged;
                    break;
                }
                else {
                    throw std::logic_error("solver_parameters.line_search_fail_action is unknown");
                }
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

        if (result->result_code != numerical_result_code_t::Converged) {
            result->score = convergence_score_t::Poor;
        }

        result->residuals_norm = residuals.objective_function(r);
        if (std::isfinite(solver_parameters.residuals_norm)) {
            if (result->residuals_norm > solver_parameters.residuals_norm) 
            {
                result->result_code = numerical_result_code_t::NotConverged;
                result->score = convergence_score_t::Poor;
            }
        }
        if (!residuals.custom_success_criteria(r, argument)) {
            result->result_code = numerical_result_code_t::CustomCriteriaFailed;
            result->score = convergence_score_t::Poor;
        }

    }

};

template <>
inline double fixed_newton_raphson<1>::argument_increment_factor(
    const double& argument,
    const double& argument_increment)
{
    double arg = std::max(1.0, abs(argument));
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
        double arg = std::max(1.0, abs(argument[component]));
        double inc = argument_increment[component];
        squared_sum += pow(inc / arg, 2);
    }
    return sqrt(squared_sum) / Dimension;
}
