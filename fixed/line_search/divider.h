#pragma once


/// @brief Настройки поиска дробления шага
struct divider_search_parameters {
    /// @brief Максимальная величина регулировки шага
    double maximum_step{ 1 };
    /// @brief Минимальная величина регулировка шага
    double minimum_step{ 0.05 };
    /// @brief Величина дробления шага
    double step_divider{ 1.5 };
    /// @brief Величина шага в случае неудачного завершения регулировки шага 
    /// Задействуется, только если это разрешено в Ньютон-Рафсоне!
    double step_on_search_fail() const {
        return minimum_step;
    }
};

/// @brief Линейный поиск методом дробления шага
class divider_search {
public:
    /// @brief Нужно для обобщения использования в Ньютоне-Рафсоне
    typedef divider_search_parameters parameters_type;
public:
    /// @brief Алгоритм регулировки шага методом дробления 
    /// @return [величина спуска; количество итераций]
    /// Величина спуска = nan, если произошел сбой регулировка
    template <typename Function>
    static inline std::pair<double, size_t> search(
        const divider_search_parameters& parameters,
        Function& function,
        double a, double b,
        double f_a, double f_b
    )
    {
        bool found_better_than_initial = false;

        double function_initial = f_a;

        double alpha_prev, alpha_curr;
        double function_prev, function_current;

        // найти максимальное значение alpha, при котором расчет не разваливается
        alpha_curr = b;// parameters.maximum_step;
        function_current = f_b;

        size_t index = 1;

        while (alpha_curr >= parameters.minimum_step)
        {
            alpha_prev = alpha_curr;
            alpha_curr /= parameters.step_divider;
            function_prev = function_current;
            function_current = function(alpha_curr);

            if (function_current < function_initial)
                found_better_than_initial = true;

            if (found_better_than_initial) {
                // ждем, когда ц.ф. начнет расти
                if (function_current > function_prev) {
                    return { alpha_prev, index };
                }
                else
                    continue;
            }
            index++;
        }

        return { std::numeric_limits<double>::quiet_NaN(), index };
    }
};

