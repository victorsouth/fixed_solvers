#pragma once
#include "golden_section.h"

/// @brief Параметры поиска золотого сечения с учётом области определения.
struct golden_section_domain_discovery_parameters : golden_section_parameters {
    /// @brief Режим обработки области определения
    domain_discovery_mode_t mode{ domain_discovery_mode_t::require_connected_domain };
    /// @brief Если true, проверка унимодальности по локальному максимуму не выполняется
    bool allow_non_unimodal{ false };
};

/// @brief Оптимизация методом золотого сечения с учётом ограниченной ООФ
class golden_section_search_domain_discovery {
public:
    /// @brief Для использования из Ньютона-Рафсона
    typedef golden_section_domain_discovery_parameters parameters_type;
public:
    /// @brief Оптимизация методом золотого сечения с обработкой выхода за область определения функции
    /// @details Выход за ООФ маркируется исключением domain_violation. Если f(b) передано конечным,
    /// оно используется без вызова в b; иначе (по умолчанию quiet_NaN) — evaluate(b).
    /// @return [величина спуска; количество итераций]
    /// величина спуска = nan, если произошел сбой регулировки
    template <typename Function>
    static inline std::pair<double, size_t> search(
        const golden_section_domain_discovery_parameters& parameters,
        Function& function,
        double a, double b,
        double f_a,
        double f_b = std::numeric_limits<double>::quiet_NaN()
    )
    {
        struct evaluation_t {
            bool in_domain{ false };
            bool has_nan{ false };
            double value{ std::numeric_limits<double>::quiet_NaN() };
        };

        auto check_convergence = [&](double f_min, double f_0) {
            return parameters.decrement_factor_criteria(f_min, f_0)
                || parameters.target_value_criteria(f_min);
        };

        bool seen_domain_gap = false;

        // Оценка функции в точке
        auto evaluate = [&](double x) {
            evaluation_t result;
            try {
                result.value = function(x);
                result.in_domain = true;
                result.has_nan = !std::isfinite(result.value);
            }
            catch (const domain_violation&) {
                if (parameters.mode == domain_discovery_mode_t::forbid_exit) {
                    throw std::runtime_error("domain violation detected in forbid_exit mode");
                }
                result.in_domain = false;
                if (parameters.mode == domain_discovery_mode_t::require_connected_domain) {
                    if (x > a) {
                        if (seen_domain_gap) {
                            throw std::logic_error("disconnected domain detected");
                        }
                        seen_domain_gap = true;
                    }
                }
            }
            return result;
        };

        auto fail_result = [&]() {
            return std::make_pair(
                std::numeric_limits<double>::quiet_NaN(),
                parameters.iteration_count + 1
            );
        };

        if (!std::isfinite(f_a)) {
            // Левая граница всегда должна быть в ООФ
            return fail_result();
        }

        evaluation_t b_eval;
        if (std::isfinite(f_b)) {
            b_eval.in_domain = true;
            b_eval.has_nan = false;
            b_eval.value = f_b;
        }
        else {
            b_eval = evaluate(b);
        }
        if (b_eval.has_nan) {
            // Случайный nan, нет выхода за ООФ.
            return fail_result();
        }

        const double f_0 = f_a;

        if (b_eval.in_domain) {
            if (auto boundary_alpha = parameters.try_resolve_step_by_target_value(f_a, b_eval.value, a, b)) {
                return { *boundary_alpha, 0 };
            }
        }

        // Текущий рекорд минимума: стартуем с a, и учитываем b только если он в области определения.
        double best_x = a;
        double best_f = f_a;
        if (b_eval.in_domain && b_eval.value < best_f) {
            best_x = b;
            best_f = b_eval.value;
        }

        if (check_convergence(best_f, f_0)) {
            return { best_x, 0 };
        }

        auto update_minimum = [&](
            double alpha, const evaluation_t& alpha_eval,
            double beta, const evaluation_t& beta_eval
        ) {
            // Обновление рекорда — только по определённым точкам
            if (alpha_eval.in_domain && (alpha_eval.value < best_f)) {
                best_x = alpha;
                best_f = alpha_eval.value;
            }
            if (beta_eval.in_domain && (beta_eval.value < best_f)) {
                best_x = beta;
                best_f = beta_eval.value;
            }
        };

        for (size_t index = 0; index < parameters.iteration_count; ++index) {
            const double alpha = golden_section_search::get_alpha(a, b);
            const double beta = golden_section_search::get_beta(a, b);

            const evaluation_t alpha_eval = evaluate(alpha);
            const evaluation_t beta_eval = evaluate(beta);

            if (alpha_eval.has_nan || beta_eval.has_nan) {
                return fail_result();
            }

            update_minimum(alpha, alpha_eval, beta, beta_eval);

            auto check_local_max = [&](bool allow_non_unimodal) {
                if (allow_non_unimodal) {
                    return;
                }
                if (!alpha_eval.in_domain || !beta_eval.in_domain) {
                    return;
                }
                // Проверка унимодальности только по доступным точкам.
                if (alpha_eval.value > f_a && alpha_eval.value > beta_eval.value) {
                    throw std::logic_error("non-unimodal function detected, f_alpha is max");
                }
                if (b_eval.in_domain && beta_eval.value > alpha_eval.value && beta_eval.value > b_eval.value) {
                    throw std::logic_error("non-unimodal function detected, f_beta is max");
                }
            };

            check_local_max(parameters.allow_non_unimodal);

            const bool defined_b = b_eval.in_domain;
            const bool defined_alpha = alpha_eval.in_domain;
            const bool defined_beta = beta_eval.in_domain;

            const bool allow_heuristic =
                parameters.mode == domain_discovery_mode_t::allow_disconnected_domain;

            if (!defined_alpha && !defined_beta) {
                // случаи 15, 16 — нет информации из унимодальности
                if (!allow_heuristic) {
                    return fail_result();
                }
                b = alpha;
                b_eval = alpha_eval;
            }
            else if (!defined_alpha) {
                // defined_beta = true, случаи 8-10
                if (defined_b && (beta_eval.value > b_eval.value)) {
                    // случай 8
                    a = beta;
                    f_a = beta_eval.value;
                }
                else if (beta_eval.value < f_a) {
                    // случай 9: новый отрезок [a, beta]
                    b = beta;
                    b_eval = beta_eval;
                }
                else {
                    // случай 10
                    throw std::logic_error("minimum is not in [a, b]");
                }
            }
            else if (!defined_beta) {
                if (!defined_b) {
                    // случаи 11-14: эвристика зависит от флага.
                    // При отключенной эвристике (случай 12) нужно запускать полное исследование ООФ.
                    // В текущем контракте метода это отражаем fail-результатом.
                    if (!allow_heuristic) {
                        return fail_result();
                    }
                    if (alpha_eval.value < f_a) {
                        // случай 13: новый отрезок [a, beta]
                        b = beta;
                        b_eval = beta_eval;
                    }
                    else {
                        // случай 11: новый отрезок [a, alpha]
                        b = alpha;
                        b_eval = alpha_eval;
                    }
                }
                else {
                    // defined_b = true, случаи 5-7
                    if (alpha_eval.value < f_a) {
                        // случай 5
                        b = alpha;
                        b_eval = alpha_eval;
                    }
                    else if (alpha_eval.value > b_eval.value) {
                        // случай 6
                        a = alpha;
                        f_a = alpha_eval.value;
                    }
                    else {
                        // случай 7
                        throw std::logic_error("minimum is not in [a, b]");
                    }
                }
            }
            else {
                // defined_alpha = true, defined_beta = true, случаи 1-4
                if (alpha_eval.value < beta_eval.value) {
                    // случаи 1, 3
                    b = beta;
                    b_eval = beta_eval;
                }
                else {
                    // случаи 2, 4
                    a = alpha;
                    f_a = alpha_eval.value;
                }
            }

            if (check_convergence(best_f, f_0)) {
                return { best_x, index + 1 };
            }
        }

        if (best_f < f_0) {
            return { best_x, parameters.iteration_count + 1 };
        }

        return fail_result();
    }
};
