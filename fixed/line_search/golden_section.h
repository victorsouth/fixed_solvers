#pragma once

/// @brief Исключение - выход за область определения функции
struct domain_violation {};

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

    /// @brief Оптимизация методом золотого сечения с обработкой выхода за область определения функции
    /// @details Выход за ООФ маркируется исключением domain_violation
    /// @return [величина спуска; количество итераций]
    /// величина спуска = nan, если произошел сбой регулировки
    template <typename Function>
    static inline std::pair<double, size_t> search_domain_based(
        const golden_section_parameters& parameters,
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

        // Оценка функции в точке
        auto evaluate = [&](double x) {
            evaluation_t result;
            try {
                result.value = function(x);
                result.in_domain = true;
                result.has_nan = !std::isfinite(result.value);
            }
            catch (const domain_violation&) {
                result.in_domain = false;
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
        auto b_eval = std::isfinite(f_b)
            ? evaluation_t{ true, false, f_b }
            : evaluate(b);
        if (b_eval.has_nan) {
            // Случайный nan, нет выхода за ООФ.
            return fail_result();
        }
        if (b_eval.in_domain) {
            f_b = b_eval.value;
        }

        double f_0 = f_a;
        auto [x_min, f_min] = (b_eval.in_domain && b_eval.value < f_a)
            ? std::make_pair(b, b_eval.value)
            : std::make_pair(a, f_a);

        if (check_convergence(f_min, f_0)) {
            return { x_min, 0 };
        }

        double alpha = get_alpha(a, b);
        double beta = get_beta(a, b);

        auto alpha_eval = evaluate(alpha);
        auto beta_eval = evaluate(beta);
        if (alpha_eval.has_nan || beta_eval.has_nan) {
            return fail_result();
        }

        auto check_local_max = [&]() {
            if (alpha_eval.in_domain && beta_eval.in_domain) {
                // Проверка унимодальности только по доступным точкам.
                if (alpha_eval.value > f_a && alpha_eval.value > beta_eval.value) {
                    throw std::logic_error("non-unimodal function detected, f_alpha is max");
                }
                if (b_eval.in_domain && beta_eval.value > alpha_eval.value && beta_eval.value > b_eval.value) {
                    throw std::logic_error("non-unimodal function detected, f_beta is max");
                }
            }
        };

        auto update_minimum = [&]() {
            const bool defined_alpha = alpha_eval.in_domain;
            const bool defined_beta = beta_eval.in_domain;

            if (!defined_alpha) {
                if (defined_beta) {
                    // разрывная ООФ. Для выбора стороны сравниваем f(beta) и f(a).
                    // Если beta не хуже a, рекорд переносим в beta, иначе оставляем у a.
                    std::tie(x_min, f_min) = (beta_eval.value <= f_a)
                        ? std::make_pair(beta, beta_eval.value)
                        : std::make_pair(a, f_a);
                }
                else {
                // alpha и beta вне ООФ. Рекорд остается у a.
                std::tie(x_min, f_min) = std::make_pair(a, f_a);
                }
            }
            else if (!defined_beta) {
                // beta вне ООФ. Сравниваем только alpha и a.
                std::tie(x_min, f_min) = (alpha_eval.value < f_a)
                    ? std::make_pair(alpha, alpha_eval.value)
                    : std::make_pair(a, f_a);
            }
            else {
                // alpha и beta в ООФ.
                if (alpha_eval.value < beta_eval.value) {
                    std::tie(x_min, f_min) = (alpha_eval.value < f_a)
                        ? std::make_pair(alpha, alpha_eval.value)
                        : std::make_pair(a, f_a);
                }
                else if (b_eval.in_domain) {
                    // Текущий рекорд либо beta, либо b
                    std::tie(x_min, f_min) = (b_eval.value < beta_eval.value)
                        ? std::make_pair(b, b_eval.value)
                        : std::make_pair(beta, beta_eval.value);
                }
                else {
                    // Рекорд остается у beta.
                    std::tie(x_min, f_min) = std::make_pair(beta, beta_eval.value);
                }
            }
        };

        check_local_max();
        update_minimum();
        if (check_convergence(f_min, f_0)) {
            return { x_min, 1 };
        }

        for (size_t index = 1; index < parameters.iteration_count; ++index) {
            const bool defined_alpha = alpha_eval.in_domain;
            const bool defined_beta = beta_eval.in_domain;

            if (!defined_alpha) {
                if (defined_beta) {
                    // Разрыв ООФ. Эвристика минимума функции
                    // 1) если f(beta) <= f(a), переходим в [beta, b];
                    // 2) иначе переходим в [a, alpha].
                    if (beta_eval.value <= f_a) {
                        a = beta;
                        f_a = beta_eval.value;
                    }
                    else {
                        b = alpha;
                        b_eval = alpha_eval;
                    }
                    alpha = get_alpha(a, b);
                    alpha_eval = evaluate(alpha);
                    beta = get_beta(a, b);
                    beta_eval = evaluate(beta);
                    continue;
                }

                // Полагаем, что [alpha, b] не в ООФ. Ищем минимум в [a, alpha]
                b = alpha;
                b_eval = alpha_eval;
                alpha = get_alpha(a, b);
                alpha_eval = evaluate(alpha);
                beta = get_beta(a, b);
                beta_eval = evaluate(beta);
            }

            else {
                if (!defined_beta) {
                    // Полагаем, что [beta, b] не в ООФ. Ищем минимум в [a, beta]
                    b = beta;
                    b_eval = beta_eval;
                    beta = alpha;
                    beta_eval = alpha_eval;
                    alpha = get_alpha(a, b);
                    alpha_eval = evaluate(alpha);
                }
                else {
                    // alpha и beta в ООФ. Стандартная итерация ЗС.
                    if (alpha_eval.value < beta_eval.value) {
                        // минимум в [a, beta]
                        b = beta;
                        if (b_eval.in_domain) {
                            f_b = b_eval.value;
                        }
                        b_eval = beta_eval;
                        beta = alpha;
                        beta_eval = alpha_eval;
                        alpha = get_alpha(a, b);
                        alpha_eval = evaluate(alpha);
                    }
                    else {
                        // минимум в [alpha, b]
                        a = alpha;
                        f_a = alpha_eval.value;
                        alpha = beta;
                        alpha_eval = beta_eval;
                        beta = get_beta(a, b);
                        beta_eval = evaluate(beta);
                    }
                }
            }

            if (alpha_eval.has_nan || beta_eval.has_nan) {
                return fail_result();
            }

            check_local_max();
            update_minimum();

            if (check_convergence(f_min, f_0)) {
                return { x_min, index + 1 };
            }
        }

        if (f_min < f_0) {
            return { x_min, parameters.iteration_count + 1 };
        }
        
        return fail_result();
    }
};

