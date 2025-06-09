/*!
* \file
* \brief В данном .h файле реализован решатель Ньютона - Рафсона
* 
* \author Автор файла - В.В. Южанин, Автор документации - И.Б. Цехместрук
* \date Дата докуменатции - 2023-03-31
* 
* Документация к этому файлу находится в:
* 1. 01. Антуан Рауль-Дальтон\02. Документы - черновики\Иван\01. Описание численного метода
* 2. 01. Антуан Рауль-Дальтон\04. Внутренние учебные материалы\Метод Ньютона-Рафсона
*/

#pragma once

#include "line_search/golden_section.h"
#include "line_search/divider.h"

#include <iostream>
#include <vector>
#include <map>
#include <iostream>
#include <string>
#include <algorithm>

// Подключаемые библиотеки
using std::vector;
using std::map;
using std::wstringstream;
using std::wstring;

enum class step_constraint_algorithm_t { Quadprog, CoordinateDescent };



///@brief Параметры диагностики сходимости
struct fixed_solver_analysis_parameters_t {
    ///@brief Собирает историю значений аргумента
    bool argument_history{ false };
    /// @brief История целевой функции
    bool objective_function_history{ false };
    ///@brief Собирает значения шага в методе Ньютона-Рафсона
    bool steps{false};
    ///@brief Выполнятся исследование целевой функции на каждом шаге
    bool line_search_explore{ false };
};

///@brief Действия при неудачной регулировке шага
enum class line_search_fail_action_t { 
    TreatAsFail, ///< Считать ошибкой 
    TreatAsSuccess, ///< Считать, что все нормально
    PerformMinStep ///< Выбрать минмимальный шаг
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
    /// @brief Использовать условную оптимизацию при необходимости ограничить шаг
    bool step_constraint_as_optimization{ false };
    /// @brief Алгоритм условной оптимизации для ограничения шага
    step_constraint_algorithm_t step_constraint_algorithm{ step_constraint_algorithm_t::CoordinateDescent };
    /// @brief Линейные ограничения
    fixed_linear_constraints<Dimension, LinearConstraintsCount> linear_constraints;
    /// @brief Параметры алгоритма регулировки шага
    typename LineSearch::parameters_type line_search;
    /// @brief Действие при ошибке регулировки шага
    line_search_fail_action_t line_search_fail_action{ line_search_fail_action_t::TreatAsFail };
    /// @brief Количество итераций
    size_t iteration_count{ 100 };
    /// @brief Критерий выхода из процедуры (погрешность метода) по приращению аргумента 
    double argument_increment_norm{ 1e-4 };
    /// @brief Критерий выхода из процедуры (погрешность метода) по норме невязок
    double residuals_norm{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief Разрешить успешный выход по невязкам, даже если алгоритм неудачно завершился 
    /// (помогает, когда уперлись в ограничение аргумента)
    bool residuals_norm_allow_force_success{ false };
    /// @brief Проверять векторный шага после регулировки или перед ней
    bool step_criteria_assuming_search_step{ false };
};

/// @brief Результат расчета Ньютона - Рафсона
enum class numerical_result_code_t
{
    NoNumericalError, IllConditionedMatrix, LargeConditionNumber, CustomCriteriaFailed,
    NotConverged, NumericalNanValues, LineSearchFailed, Converged, InProgress
};

/// \brief в поток
/// \param os - поток
/// \param c - numerical_result_code_t
/// \return поток
inline std::ostream& operator<<(std::ostream& os, const numerical_result_code_t& c) {
    switch (c) {
    case numerical_result_code_t::NoNumericalError:
        os << "NoNumericalError(" << int(c) << ")";
        break;
    case numerical_result_code_t::NumericalNanValues:
        os << "NumericalNanValues(" << int(c) << ")";
        break;
    case numerical_result_code_t::Converged:
        os << "Converged(" << int(c) << ")";
        break;
    case numerical_result_code_t::NotConverged:
        os << "NotConverged(" << int(c) << ")";
        break;
    case numerical_result_code_t::CustomCriteriaFailed:
        os << "CustomCriteriaFailed(" << int(c) << ")";
        break;
    case numerical_result_code_t::IllConditionedMatrix:
        os << "IllConditionedMatrix(" << int(c) << ")";
        break;
    case numerical_result_code_t::LargeConditionNumber:
        os << "LargeConditionNumber(" << int(c) << ")";
        break;
    case numerical_result_code_t::LineSearchFailed:
        os << "LineSearchFailed(" << int(c) << ")";
        break;
    }
    return os;
}


/// @brief Оценка качества сходимости
enum class convergence_score_t : int {
    Excellent = 5, Good = 4, Satisfactory = 3, Poor = 2, Error = 1
};

/// @brief Преобразование кода сходимости в строковый вид (wide string)
inline const map<convergence_score_t, wstring>& get_score_wstrings() {
    static const map<convergence_score_t, wstring> score_strings{
        {convergence_score_t::Excellent, L"Excellent(5)"},
        {convergence_score_t::Good , L"Good(4)"},
        {convergence_score_t::Satisfactory, L"Satisfactory(3)"},
        {convergence_score_t::Poor, L"Poor(2)"},
        {convergence_score_t::Error, L"Error(1)"},
    };
    return score_strings;
}

/// @brief Преобразование кода сходимости в строковый вид (regular string)
inline const map<convergence_score_t, std::string>& get_score_strings() {
    static const map<convergence_score_t, std::string> score_strings{
        {convergence_score_t::Excellent, "Excellent(5)"},
        {convergence_score_t::Good , "Good(4)"},
        {convergence_score_t::Satisfactory, "Satisfactory(3)"},
        {convergence_score_t::Poor, "Poor(2)"},
        {convergence_score_t::Error, "Error(1)"},
    };
    return score_strings;
}



/// \brief в поток
/// \param os - поток
/// \param c - convergence_score_t
/// \return поток
inline std::ostream& operator<<(std::ostream& os,const convergence_score_t& c){
    os << get_score_strings().at(c);
    return os;
}

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
    const auto& score_strings = get_score_wstrings();
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
    /// @brief Тип аргумента
    typedef typename fixed_system_types<Dimension>::var_type var_type;
    /// @brief Тип функции
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
template <std::ptrdiff_t Dimension>
struct fixed_solver_result_analysis_t {
public:
    /// @brief Значения целевой функции для одной регулировки шага
    typedef vector<double> target_function_values_t;
    /// @brief Тип аргумента
    typedef typename fixed_system_types<Dimension>::var_type var_type;
    /// @brief Тип функции
    typedef typename fixed_system_types<Dimension>::right_party_type function_type;
    /// @brief Тип коэффициентов уравнения (неясно, зачем нужно)
    typedef typename fixed_system_types<Dimension>::equation_coeffs_type equation_coeffs_type;
public:
    /// @brief Исследование целевой функции по всем шагам Ньютона-Рафсона или Гаусса-Ньютона
    vector<target_function_values_t> target_function;
    /// @brief Результат расчета
    vector<var_type> argument_history;
    /// @brief Величина шагов Ньютона-Рафсона
    vector<double> steps;
    /// @brief Строит результат на основе собранного в target_function в режиме однократного расчета ц.ф.
    vector<double> get_learning_curve() const {
        vector<double> result;
        result.reserve(target_function.size());
        std::transform(target_function.begin(), target_function.end(), 
            std::back_inserter(result),
            [](const vector<double>& objective_function_value) {
                if (objective_function_value.size() != 1)
                    throw std::runtime_error("Unexpected objective function values per step");
                return objective_function_value[0];
            }
            );

        return result;
    }
};

/// @brief Настройки регулировки шага "без регулировки"
struct no_line_search_parameters {
    /// @brief Максимальная величина шага (формально нужна)
    double maximum_step{ 1.0 };
    /// @brief Ошибка при вылете регулировки 
    /// (по идее, не должна вызываться, т.к. данная "регулировка" никогда не фейлится)
    double step_on_search_fail() const {
        return 1.0;
    }
};

/// @brief Алгоритм без регулировки шага
class no_line_search {
public:
    /// @brief Формально надо для вызова из Ньютона-Рафсона
    typedef no_line_search_parameters parameters_type;

public:
    /// @brief Фиктивная регулировка шага
    /// @return [величина спуска == 1.0; количество итераций == 1]
    template <typename Function>
    static inline std::pair<double, size_t> search(
        const no_line_search_parameters& parameters,
        Function& function,
        double a, double b,
        double f_a, double f_b
    )
    {
        size_t index = 1;
        double step = 1;
        return { step, index };
    }
};




/// @brief Солвер Ньютона-Рафсона для систем уравнений фиксированной размерности
/// В том числе скалярных
/// Ньютона-Рафсона - означает регулировку шага
/// Реализована регулировка шага за счет дробления
template <std::ptrdiff_t Dimension>
class fixed_newton_raphson {
public:
    /// @brief Тип аргумента
    typedef typename fixed_system_types<Dimension>::var_type var_type;
    /// @brief Тип функции
    typedef typename fixed_system_types<Dimension>::right_party_type function_type;
    /// @brief Тип коэффициентов уравнения (неясно, зачем нужно)
    typedef typename fixed_system_types<Dimension>::equation_coeffs_type equation_coeffs_type;
private:

    /// @brief Проверка значения на Nan/infinite для скалярного случая
    /// @param value Проверяемое значение
    /// @return true, если найдены nan
    static inline bool has_not_finite(const double value)
    {
        if (!std::isfinite(value)) {
            return true;
        }
        return false;
    }

    /// @brief Поиск nan-значений в векторе
    /// @param values Проверяемое значение 
    /// @return true, если найдены nan
    static inline bool has_not_finite(const VectorXd& values)
    {
        for (ptrdiff_t index = 0; index < values.size(); ++index) {
            if (!std::isfinite(values(index))) {
                return true;
            }
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
    /// @brief Расчет относительного приращения (реализация зависит от var_type)
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

        return target_function;
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

    template <std::ptrdiff_t LinearConstraintsCount, typename LineSearch>
    static double solve_coodinate_descent(
        const fixed_solver_parameters_t<Dimension, LinearConstraintsCount, LineSearch>& solver_parameters,
        fixed_system_t<Dimension>& residuals,
        const var_type& current_residuals_value, const var_type& argument,
        size_t var_index
    )
    {
        if constexpr (Dimension == -1) {
            vector<Eigen::Triplet<double>> Jcol = residuals.jacobian_sparse_column(argument, var_index);

            SparseMatrix<double> J(current_residuals_value.size(), 1);
            J.setFromTriplets(Jcol.begin(), Jcol.end());

            Eigen::JacobiSVD<MatrixXd> svd_solver(
                J, Eigen::ComputeThinU | Eigen::ComputeThinV);
            VectorXd var_result = svd_solver.solve(-current_residuals_value);
            return var_result(0);
        }
        else if constexpr (Dimension > 1) {
            function_type Jcol = residuals.jacobian_column(argument, var_index);
            MatrixXd J = MatrixXd::Map(&Jcol[0], Dimension, 1);
            Eigen::JacobiSVD<MatrixXd> svd_solver(
                J, Eigen::ComputeThinU | Eigen::ComputeThinV);

            VectorXd var_result = svd_solver.solve(-VectorXd::Map(&current_residuals_value[0], Dimension));
            return var_result(0);
        }
        else {
            throw std::runtime_error("Coordinate descent for dimension = 1 is senseless");
        }


    }
    template <std::ptrdiff_t LinearConstraintsCount, typename LineSearch>
    static var_type solve_quadprog(
        const fixed_solver_parameters_t<Dimension, LinearConstraintsCount, LineSearch>& solver_parameters,
        fixed_system_t<Dimension>& residuals,
        const var_type& current_residuals_value, const var_type& argument
    )
    {

        if constexpr (Dimension == -1) {
#ifndef FIXED_USE_QP_SOLVER
            throw std::runtime_error("Compiled without QP support");
#endif
#ifdef FIXED_USE_QP_SOLVER
            vector<Eigen::Triplet<double>> J_triplets = residuals.jacobian_sparse(argument);

            SparseMatrix<double> J(argument.size(), argument.size());
            J.setFromTriplets(J_triplets.begin(), J_triplets.end());

            SparseMatrix<double> H = J.transpose() * J;
            VectorXd f = current_residuals_value.transpose() * J; // r * J

            auto [min, max] = solver_parameters.constraints.get_relative_constraints(argument);

            VectorXd result = solve_quadprog_box(H, f, min, max);
            return result;
#endif
        }
        else if constexpr (Dimension != 1) {
#ifndef FIXED_USE_QP_SOLVER
            throw std::runtime_error("Compiled without QP support");
#endif
#ifdef FIXED_USE_QP_SOLVER

            // Расчет Якобиана
            auto J = residuals.jacobian_sparse(argument);

            SparseMatrix<double> Js(Dimension, Dimension);
            Js.setFromTriplets(J.begin(), J.end());

            SparseMatrix<double> H = Js.transpose() * Js;
            VectorXd r = VectorXd::Map(current_residuals_value.data(), current_residuals_value.size());
            VectorXd f = r.transpose() * Js; // r * J

            //SparseMatrix<double> A;
            //VectorXd b;
            //VectorXd arg = VectorXd::Map(argument.data(), argument.size());
            //std::tie(A, b) = solver_parameters.constraints.get_inequalities_constraints_sparse(argument);
            auto [min, max] = solver_parameters.constraints.get_relative_constraints(argument);

            VectorXd result = solve_quadprog_box(H, f, min, max);

            var_type p;
            std::copy(&result(0), &result(0) + result.size(), p.begin());
            return p;
#endif
        }
        else
        {
            throw std::runtime_error("Dimension=1 should not be called with quadprog");
        }

    }

    static auto solve_newton(
        fixed_system_t<Dimension>& residuals,
        const var_type& current_residuals_value, const var_type& argument)
    {
        if constexpr (Dimension == -1) {
            auto J_triplets = residuals.jacobian_sparse(argument);

            SparseMatrix<double> J(argument.size(), argument.size());
            J.setFromTriplets(J_triplets.begin(), J_triplets.end());

            Eigen::SparseLU<SparseMatrix<double> > solver;
            solver.analyzePattern(J);
            solver.factorize(J);
            if (solver.info() == Eigen::Success) {
                VectorXd result = -solver.solve(current_residuals_value);
                return result;
            }
            else {
                throw std::runtime_error("cannot calc newton step");
            }
        }
        else {
            // Расчет Якобиана
            auto J = residuals.jacobian_dense(argument);
            // Расчет текущего шага Ньютона
            auto result = -solve_linear_system(J, current_residuals_value);
            return result;
        }
    }


    /// @brief 
    /// @param solver_parameters 
    /// @param residuals 
    /// @param result 
    /// @param analysis 
    /// @return Истина, если расчет надо завершать 
    /// Успешно или наоборот, ошибка расчета, но завершать
    template <std::ptrdiff_t LinearConstraintsCount, typename LineSearch>
    static bool perform_vector_step(
        const fixed_solver_parameters_t<Dimension, LinearConstraintsCount, LineSearch>& solver_parameters,
        bool optimization_step,
        fixed_system_t<Dimension>& residuals,
        fixed_solver_result_t<Dimension>* result,
        fixed_solver_result_analysis_t<Dimension>* analysis
    )
    {
        var_type& r = result->residuals;
        function_type& argument = result->argument;
        double& argument_increment_metric = result->argument_increment_metric;
        bool& argument_increment_criteria = result->argument_increment_criteria;

        var_type p = optimization_step
            ? solve_quadprog(solver_parameters, residuals, r, argument)
            : solve_newton(residuals, r, argument);

        // Проверка критерия выхода по малому относительному приращению
        if (solver_parameters.step_criteria_assuming_search_step == false)
        {
            argument_increment_metric = argument_increment_factor(argument, p);
            argument_increment_criteria =
                argument_increment_metric < solver_parameters.argument_increment_norm;
            if (argument_increment_criteria) {
                result->result_code = numerical_result_code_t::Converged;
                return true;
            }
        }

        // Корректировка шага в соответствии с ограничениями
        solver_parameters.linear_constraints.trim(argument, p);
        solver_parameters.constraints.trim_max(argument, p);
        solver_parameters.constraints.trim_min(argument, p);
        solver_parameters.constraints.trim_relative(p);

        // Информация о том, как изменялась целевая функция по траектории шага
        if (analysis != nullptr && solver_parameters.analysis.line_search_explore) {
            analysis->target_function.push_back(perform_step_research(residuals, argument, p));
        }

        // Расчет корректироки шага с помощью Рафсона
        double search_step = perform_line_search<LineSearch>(
            solver_parameters.line_search, residuals, argument, r, p);
        if (analysis != nullptr && solver_parameters.analysis.steps) {
            analysis->steps.push_back(search_step);
        }

        if (std::isfinite(search_step)) {
            if (search_step < small_step_threshold) {
                // снижаем до 4-х баллов, если он выше
                result->score = std::min(result->score, convergence_score_t::Good);
            }
        }
        else {
            if (solver_parameters.line_search_fail_action == line_search_fail_action_t::PerformMinStep) {
                // понижаем балл до удовлетворительно, считаем дальше
                result->score = std::min(result->score, convergence_score_t::Satisfactory);
                search_step = solver_parameters.line_search.step_on_search_fail();
            }
            else if (solver_parameters.line_search_fail_action == line_search_fail_action_t::TreatAsFail) {
                result->result_code = numerical_result_code_t::LineSearchFailed;
                return true;
            }
            else if (solver_parameters.line_search_fail_action == line_search_fail_action_t::TreatAsSuccess) {
                result->result_code = numerical_result_code_t::Converged;
                return true;
            }
            else {
                throw std::logic_error("solver_parameters.line_search_fail_action is unknown");
            }
        }

        // Корректировка шага с помощью Рафсона
        var_type argument_increment = search_step * p;
        argument += argument_increment;

        if (analysis != nullptr && solver_parameters.analysis.argument_history) {
            analysis->argument_history.push_back(argument);
        }

        // Расчет невязок
        r = residuals.residuals(argument);
        result->residuals_norm = residuals.objective_function(r);
        if (has_not_finite(r)) {
            r = residuals.residuals(argument); // для отладки
            result->result_code = numerical_result_code_t::NumericalNanValues;
            return true;
        }

        bool custom_criteria = residuals.custom_success_criteria(r, argument);
        if (custom_criteria) {
            result->result_code = numerical_result_code_t::Converged;
            return true;
        }

        // Проверка критерия выхода по малому относительному приращению
        argument_increment_metric = solver_parameters.step_criteria_assuming_search_step
            ? argument_increment_factor(argument, argument_increment)
            : argument_increment_factor(argument, p);

        argument_increment_criteria =
            argument_increment_metric < solver_parameters.argument_increment_norm;

        result->result_code = argument_increment_criteria
            ? numerical_result_code_t::Converged
            : numerical_result_code_t::InProgress;

        return argument_increment_criteria;
    }

    /// @brief 
    /// @tparam LineSearch 
    /// @tparam LinearConstraintsCount 
    /// @param solver_parameters 
    /// @param residuals 
    /// @param result 
    /// @param analysis 
    /// @return Истина, если расчет надо завершать 
    /// Успешно или наоборот, ошибка расчета, но завершать
    template <std::ptrdiff_t LinearConstraintsCount, typename LineSearch>
    static bool perform_coodinate_descent_step(
        const fixed_solver_parameters_t<Dimension, LinearConstraintsCount, LineSearch>& solver_parameters,
        fixed_system_t<Dimension>& residuals,
        fixed_solver_result_t<Dimension>* result,
        fixed_solver_result_analysis_t<Dimension>* analysis
    )
    {

        bool has_succeded_search_step = false;

        var_type& r = result->residuals;
        function_type& argument = result->argument;
        double& argument_increment_metric = result->argument_increment_metric;
        bool& argument_increment_criteria = result->argument_increment_criteria;

        var_type p;
        if constexpr (Dimension == -1) {
            p = VectorXd::Zero(argument.size());
        }
        else if constexpr (Dimension > 1) {
            std::for_each(p.begin(), p.end(), [](double& value) { value = 0; });
        }

        argument_increment_metric = 0; // обнуляем, т.к. в этой переменной накапливаем максимум по приращениям всех переменных
        size_t substep_count = argument.size();

        size_t substeps_performed = 0;

        for (size_t substep = 0; substep < substep_count; ++substep) {
            if (substep != 0) {
                p[substep - 1] = 0.0;
            }

            double& var_p = p[substep];
            var_p = solve_coodinate_descent(solver_parameters, residuals, r, argument, substep);

            // Корректировка шага в соответствии с ограничениями
            solver_parameters.linear_constraints.trim(argument, p);
            solver_parameters.constraints.trim_max(argument, p);
            solver_parameters.constraints.trim_min(argument, p);
            solver_parameters.constraints.trim_relative(p);

            //double var_increment_metric = argument_increment_factor(argument, p);
            argument_increment_metric = std::max(abs(var_p), argument_increment_metric);

            // Информация о том, как изменялась целевая функция по траектории шага
            if (analysis != nullptr && solver_parameters.analysis.line_search_explore) {
                analysis->target_function.push_back(perform_step_research(residuals, argument, p));
            }

            // Расчет корректироки шага с помощью Рафсона
            double search_step = perform_line_search<LineSearch>(
                solver_parameters.line_search, residuals, argument, r, p);
            if (analysis != nullptr && solver_parameters.analysis.steps) {
                analysis->steps.push_back(search_step);
            }

            if (std::isfinite(search_step)) {
                if (search_step < small_step_threshold) {
                    // снижаем до 4-х баллов, если он выше
                    result->score = std::min(result->score, convergence_score_t::Good);
                }
                has_succeded_search_step = true;
            }
            else {
                if (solver_parameters.line_search_fail_action == line_search_fail_action_t::PerformMinStep) {
                    // понижаем балл до удовлетворительно, считаем дальше
                    result->score = std::min(result->score, convergence_score_t::Satisfactory);
                    search_step = solver_parameters.line_search.step_on_search_fail();
                    has_succeded_search_step = true;
                }
                else {
                    continue;
                }
            }

            ++substeps_performed;

            // Корректировка шага с помощью Рафсона
            var_type argument_increment = search_step * p;
            argument += argument_increment;


            if (analysis != nullptr && solver_parameters.analysis.argument_history) {
                analysis->argument_history.push_back(argument);
            }

            // Расчет невязок
            r = residuals.residuals(argument);
            if (has_not_finite(r)) {
                r = residuals.residuals(argument); // для отладки
                result->result_code = numerical_result_code_t::NumericalNanValues;
                return true;
            }

            bool custom_criteria = residuals.custom_success_criteria(r, argument);
            if (custom_criteria) {
                result->result_code = numerical_result_code_t::Converged;
                return true;
            }

            // Проверка критерия выхода по малому относительному приращению
            argument_increment_metric = std::max(abs(var_p), argument_increment_metric);

            double var_increment_metric = solver_parameters.step_criteria_assuming_search_step
                ? search_step * abs(var_p)
                : abs(var_p);

            argument_increment_metric = std::max(var_increment_metric, argument_increment_metric);
        }

        result->residuals_norm = residuals.objective_function(r);

        argument_increment_criteria =
            argument_increment_metric < solver_parameters.argument_increment_norm;
        if (argument_increment_criteria) {
            result->result_code = numerical_result_code_t::Converged;
            return true;
        }

        if (has_succeded_search_step == false) {
            switch (solver_parameters.line_search_fail_action)
            {
            case line_search_fail_action_t::TreatAsFail:
                result->result_code = numerical_result_code_t::LineSearchFailed;
                return true;
            case line_search_fail_action_t::TreatAsSuccess:
                result->result_code = numerical_result_code_t::Converged;
                return true;
            default:
                // Остается опция "PerformMinStep", но тогда has_succeded_search_step = true
                throw std::runtime_error("Wrong branch");
            }
        }
        else {
            result->result_code = numerical_result_code_t::InProgress;
            return false;
        }

    }

public:
    /// @brief Запуск численного метода
    /// @tparam LineSearch Алгоритм регулировки шага поиска
    /// @param residuals Функция невязок
    /// @param initial_argument Начальное приближение
    /// @param solver_parameters Настройки поиска
    /// @param result Результаты расчета
    template <
        std::ptrdiff_t LinearConstraintsCount = 0,
        typename LineSearch = divider_search>
    static void solve(
        fixed_system_t<Dimension>& residuals,
        const var_type& initial_argument,
        const fixed_solver_parameters_t<Dimension, LinearConstraintsCount, LineSearch>& solver_parameters,
        fixed_solver_result_t<Dimension>* result,
        fixed_solver_result_analysis_t<Dimension>* analysis = nullptr
    )
    {
        using std::max;
        using std::min;

        var_type& r = result->residuals;
        function_type& argument = result->argument;
        double& argument_increment_metric = result->argument_increment_metric;
        bool& argument_increment_criteria = result->argument_increment_criteria;


        argument = initial_argument;
        if (analysis != nullptr && solver_parameters.analysis.argument_history) {
            analysis->argument_history.push_back(argument);
        }

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

        // Запуск расчета методом Ньютона - Рафсона
        for (iteration = 0; iteration < solver_parameters.iteration_count; ++iteration)
        {

            bool stop_iterations;
            stop_iterations = perform_vector_step(solver_parameters, false, residuals, result, analysis);

            bool optimization_step = solver_parameters.step_constraint_as_optimization
                && solver_parameters.constraints.has_active_constraints(argument)
                && result->result_code == numerical_result_code_t::LineSearchFailed;
            if (optimization_step) {
                if constexpr (Dimension == 1) {
                    stop_iterations = perform_vector_step(solver_parameters, true, residuals, result, analysis);
                }
                else {
                    if (solver_parameters.step_constraint_algorithm == step_constraint_algorithm_t::CoordinateDescent) {
                        stop_iterations = perform_coodinate_descent_step(solver_parameters, residuals, result, analysis);
                    }
                    else {
                        stop_iterations = perform_vector_step(solver_parameters, true, residuals, result, analysis);
                    }
                }
            }

            if (stop_iterations)
                break;
        }

        // Выставление оценки от количества итераций
        if (iteration > 0.3 * solver_parameters.iteration_count) {
            result->score = std::min(result->score, convergence_score_t::Satisfactory);
        }
        else if (iteration > 0.15 * solver_parameters.iteration_count)
        {
            result->score = std::min(result->score, convergence_score_t::Good);
        }

        if (result->result_code != numerical_result_code_t::Converged) {
            result->score = convergence_score_t::Poor;
        }

        if (std::isfinite(solver_parameters.residuals_norm)) {
            if (result->residuals_norm > solver_parameters.residuals_norm)
            {
                result->result_code = numerical_result_code_t::NotConverged;
                result->score = convergence_score_t::Poor;
            }
            else {
                if (solver_parameters.residuals_norm_allow_force_success) {
                    result->result_code = numerical_result_code_t::Converged;
                    int score = max(
                        static_cast<int>(result->score),
                        static_cast<int>(convergence_score_t::Satisfactory));
                    result->score = static_cast<convergence_score_t>(score);
                }
            }
        }
    }

    /// @brief Запуск численного метода, вызывает solve
    /// это просто вызов для обратной совместимости, для тех, кто привык использовать solve_dense
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
        // это просто вызов для обратной совместимости, для тех, кто привык использовать solve_dense
        solve(residuals, initial_argument, solver_parameters, result, analysis);
    }
};

/*! @brief Расчет относительного приращения для скалярного случая
* @param argument Текущее значение
* @param argument_increment Текущее приращение
* @return Отношение приращения к текущему значению
*/
template <>
inline double fixed_newton_raphson<1>::argument_increment_factor(
    const double& argument,
    const double& argument_increment)
{
    using std::max;
    double arg = max(1.0, abs(argument));
    double inc = argument_increment;

    double result = abs(inc / arg);
    return result;
}

template <std::ptrdiff_t Dimension>
/// @brief Расчет относительного приращения для векторного случая
/// @param argument Текущее значение
/// @param argument_increment Текущее приращение
/// @return Отношение приращения к текущему значению
inline double fixed_newton_raphson<Dimension>::argument_increment_factor(
    const var_type& argument, const var_type& argument_increment)
{
    using std::max;
    double squared_sum = 0;

    for (int component = 0; component < argument_increment.size(); ++component) {
        double arg = max(1.0, abs(argument[component]));
        double inc = argument_increment[component];
        squared_sum += pow(inc / arg, 2);
    }
    if constexpr (Dimension == -1) {
        return sqrt(squared_sum) / argument.size();
    }
    else {
        return sqrt(squared_sum) / Dimension;
    }
    
}



