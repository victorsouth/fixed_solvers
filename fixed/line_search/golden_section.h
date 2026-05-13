#pragma once

/// @brief Исключение - выход за область определения функции
struct domain_violation {};

/// @brief Режимы обработки области определения при поиске золотого сечения
enum class domain_discovery_mode_t {
    /// @brief Не разрешать выход за область определения
    forbid_exit,
    /// @brief Требовать связную область определения
    require_connected_domain,
    /// @brief Разрешать несвязную область определения (использовать эвристику ООФ)
    allow_disconnected_domain
};

/// @brief Настройки поиска золотого сечения
struct golden_section_parameters {
    /// @brief Максимальный шаг расчета
    double maximum_step{ 1 };
    /// @brief Количество итераций
    size_t iteration_count{ 10 };
    /// @brief Отношение начального значения ц.ф. к итоговому
    double function_decrement_factor{ 2 };
    /// @brief Достаточное (достаточно малое) значение целевое значение ц.ф.
    /// Если NaN, то не используется
    double function_target_value{ std::numeric_limits<double>::quiet_NaN() };    
    /// @brief Если регулировка шага отвалилась и разрешено шаг все-таки сделать, то он будет такой
    double fail_step_size{ 0.05 };
    /// @brief Величина шага в случае неудачного завершения регулировки шага 
    /// Задействуется, только если это разрешено в Ньютон-Рафсоне!
    double step_on_search_fail() const {
        return fail_step_size;
    }

    /// @brief Возвращает относительную длины отрезка локализации минимума
    double get_final_section_length() const {
        return pow(0.618, iteration_count);
    }
    /// @brief Проверяет сильную сходимость по снижению ц.ф.
    bool decrement_factor_criteria(double f_current_min, double f_0) const {
        if (std::isfinite(function_decrement_factor)) {
            double decrement = f_0 / f_current_min;
            // достигнута сильная сходимость
            return  decrement > function_decrement_factor;
        }
        else
            return false;
    }
    /// @brief Проверяет достигнутое значение ц.ф.
    bool target_value_criteria(double f_current_min) const {
        return std::isfinite(function_target_value) && (f_current_min < function_target_value);
    }
};


/// @brief Оптимизация методом золотого сечения
class golden_section_search {
public:
    /// @brief Для использования из Ньютона-Рафсона
    typedef golden_section_parameters parameters_type;
public:
    /// @brief Левая граница ЗС
    static inline double get_alpha(double a, double b) {
        return a + 2 / (3 + sqrt(5)) * (b - a);
    }
    /// @brief Правая граница ЗС
    static inline double get_beta(double a, double b) {
        return a + 2 / (1 + sqrt(5)) * (b - a);
    }
public:
    /// @brief Оптимизация методом золотого сечения
    /// @return [величина спуска; количество итераций]
    /// величина спуска = nan, если произошел сбой регулировки
    template <typename Function>
    static inline std::pair<double, size_t> search(
        const golden_section_parameters& parameters,
        Function& function,
        double a, double b,
        double f_a, double f_b
    )
    {
        auto check_convergence = [&](double f_min, double f_0) {
            return parameters.decrement_factor_criteria(f_min, f_0)
                || parameters.target_value_criteria(f_min);
        };

        // Начальное значение ц.ф., по которому смотрим фактор снижения
        double f_0 = f_a;

        // Инициализируем текущий рекорд
        auto [x_min, f_min] = f_b < f_a
            ? std::make_pair(b, f_b)
            : std::make_pair(a, f_a);
        // Если снижение сразу достаточно, выходим
        if (check_convergence(f_min, f_0)) {
            // трактуем, что вообще не понадобилось итераций
            return { x_min, 0 };
        }

        // Нулевая итерация, надо рассчитать два значения f(alpha), f(beta)
        double alpha = get_alpha(a, b);
        double f_alpha = function(alpha);

        double beta = get_beta(a, b);
        double f_beta = function(beta);

        // функция проверят наличие экстремума-максимума в четрыем известных точках: {a, alpha, beta, b}
        // что противоречит унимодальности функции
        auto check_local_max = [&]() {
            if (f_alpha > f_a && f_alpha > f_beta) {
                throw std::logic_error("non-unimodal function detected, f_alpha is max");
            }
            if (f_beta > f_alpha && f_beta > f_b) {
                throw std::logic_error("non-unimodal function detected, f_beta is max");
            }
        };
        //check_local_max(); // отключим, поскольку шум

        // функция обновляет рекорд x_min, f_min по известным значениям a, alpha, beta, b
        auto update_minimum = [&x_min = x_min, &f_min = f_min, &f_alpha, &f_beta, &a, &b, &f_a, &f_b, &alpha, &beta]() {            // Выбор текущего рекорда по минимуму
            if (f_alpha < f_beta) {
                // 1. исходя из унимодальности, минимум лежит в диапазоне [a, beta]
                // 2. из этого диапазона рассчитаны значения f(a), f(alpha), f(beta)
                // нужно найти минимум из них
                // 3. при этом уже проверили, что f(alpha) < f(beta)
                // остается найти min(f(alpha), f(a))
                std::tie(x_min, f_min) = f_alpha < f_a
                    ? std::make_pair(alpha, f_alpha)
                    : std::make_pair(a, f_a);
            }
            else {
                // 1. минимум в диапазоне [alpha, b]
                // 2. из этого диапазона рассчитаны значения f(alpha), f(beta), f(b)
                // 3. при этом уже проверили, что f(alpha) > f(beta)
                // остается найти min(f(b), f(beta))
                std::tie(x_min, f_min) = f_b < f_beta
                    ? std::make_pair(b, f_b)
                    : std::make_pair(beta, f_beta);
            }
        };
        update_minimum();
        if (check_convergence(f_min, f_0)) {
            return { x_min, 1 };
        }

        // Цикл на n-1 итерацию, т.к. нулевая итерация отличается от остальных и уже сделана
        for (size_t index = 1; index < parameters.iteration_count; ++index) {
            if (f_alpha < f_beta) {
                // Переходим в [a, beta]
                // Отрезок [a, alpha] - "длинный" отрезок золотого сечения
                b = beta; // b поменяли, a остается
                f_b = f_beta;
                beta = alpha;
                f_beta = f_alpha;
                alpha = get_alpha(a, b);
                f_alpha = function(alpha);
            }
            else {
                // Переходим в [alpha, b]. 
                // Отрезок [beta, b] - "длинный" отрезок золотого сечения
                a = alpha; // a поменяли, b остается
                f_a = f_alpha;
                alpha = beta;
                f_alpha = f_beta;
                beta = get_beta(a, b);
                f_beta = function(beta);
            }
            update_minimum();
            if (check_convergence(f_min, f_0)) {
                return { x_min, index + 1 };
            }
        }

        if (f_min < f_0) {
            // не сошлись на последней итерации, но слабая сходимость есть
            // чтобы отличалось от сильной сходимости на последней итерации, делаем n + 1
            return { x_min, parameters.iteration_count + 1 };
        }
        else {
            return {
                std::numeric_limits<double>::quiet_NaN(),
                parameters.iteration_count + 1
            };
        }
    }

};

/// @brief Параметры поиска золотого сечения с учётом области определения.
struct golden_section_domain_discovery_parameters : golden_section_parameters {
    /// @brief Режим обработки области определения
    domain_discovery_mode_t mode{ domain_discovery_mode_t::require_connected_domain };
};

class golden_section_search_domain_discovery {
public:
    /// @brief Для использования из Ньютона-Рафсона
    typedef golden_section_domain_discovery_parameters parameters_type;
public:
    /// @brief Оптимизация методом золотого сечения с обработкой выхода за область определения функции
    /// @details Выход за ООФ маркируется исключением domain_violation
    /// @return [величина спуска; количество итераций]
    /// величина спуска = nan, если произошел сбой регулировки
    template <typename Function>
    static inline std::pair<double, size_t> search(
        const golden_section_domain_discovery_parameters& parameters,
        Function& function,
        double a, double b,
        double f_a, double f_b
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

        // Если f_b известна, переиспользуем ее. Иначе пробуем вычислить.
        evaluation_t b_eval = std::isfinite(f_b)
            ? evaluation_t{ true, false, f_b }
            : evaluate(b);
        if (b_eval.has_nan) {
            // Случайный nan, нет выхода за ООФ.
            return fail_result();
        }
        if (b_eval.in_domain) {
            f_b = b_eval.value;
        }

        const double f_0 = f_a;

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

            auto check_local_max = [&]() {
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

            check_local_max();

            const bool defined_b = b_eval.in_domain;
            const bool defined_alpha = alpha_eval.in_domain;
            const bool defined_beta = beta_eval.in_domain;

            // Выбор нового интервала согласно псевдокоду.
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

