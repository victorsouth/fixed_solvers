#pragma once

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

