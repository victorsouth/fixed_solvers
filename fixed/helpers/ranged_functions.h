#pragma once

#include <vector>
#include <algorithm>
#include <stdexcept>
#include <utility>
#include <cstddef>

namespace fixed_solvers {

/// @brief Проверяет, находится ли значение в заданном диапазоне.
/// @param value Проверяемое значение.
/// @param range_begin Начало диапазона.
/// @param range_end Конец диапазона.
/// @return true, если значение находится в диапазоне [range_begin, range_end] (или [range_end, range_begin], если range_begin > range_end).
inline bool value_in_range(double value, double range_begin, double range_end)
{
    if (range_begin > range_end)
        return value_in_range(value, range_end, range_begin);

    return value >= range_begin && value <= range_end;
}

/// @brief Вычисляет коэффициенты интеграла полинома.
/// @details НЕ используется формат Matlab! В отличие от функции,
/// индекс коэффициента является показателем степени.
/// @tparam Coefficients Тип контейнера коэффициентов полинома.
/// @param poly_coeffs Коэффициенты полинома.
/// @return Вектор коэффициентов интеграла полинома.
template <typename Coefficients>
inline std::vector<double> poly_integral_coefficients(const Coefficients& poly_coeffs)
{
    // НЕ используется формат Matlab! в отличие от функции 
    // индекс коэффициента является показателем степени

    const size_t n = poly_coeffs.size();
    std::vector<double> result(n + 1);

    for (size_t index = 1; index < n + 1; ++index)
    {
        result[index] = poly_coeffs[index - 1] / index;
    }
    return result;
}

/// @brief Структура, описывающая диапазон функции и её коэффициенты.
/// @details Используется для задания кусочно-определённых функций, где каждый диапазон 
/// имеет свои коэффициенты.
/// @tparam Coeffs Тип коэффициентов, используемых для аппроксимации функции.
template <typename Coeffs>
struct function_range_t {
    /// @brief Левая граница диапазона.
    double range_start;
    /// @brief Правая граница диапазона.
    double range_end;
    /// @brief Коэффициенты, определяющие поведение функции на данном диапазоне.
    Coeffs coefficients;
};

/// @brief Базовый класс для работы с функциями, заданными на множестве диапазонов.
/// @details Хранит набор диапазонов с соответствующими коэффициентами и предоставляет 
/// метод для поиска диапазона по значению аргумента.
/// @tparam Coeffs Тип коэффициентов, используемых в диапазонах.
template <typename Coeffs>
class ranged_function_t {
protected:
    /// @brief Вектор диапазонов функции.
    std::vector<function_range_t<Coeffs>> ranges;
public:
    /// @brief Конструктор по умолчанию.
    ranged_function_t() = default;
    /// @brief Оператор присваивания по умолчанию.
    ranged_function_t& operator= (const ranged_function_t&) = default;

    /// @brief Конструктор с инициализацией диапазонов.
    /// @param _ranges Вектор диапазонов функции.
    /// @details Диапазоны автоматически сортируются по начальной границе.
    ranged_function_t(const std::vector<function_range_t<Coeffs>>& _ranges)
        : ranges(_ranges)
    {
        std::sort(ranges.begin(), ranges.end(),
            [](const function_range_t<Coeffs>& r1, const function_range_t<Coeffs>& r2)
            {
                return r1.range_start < r2.range_start;
            });
    }

    /// @brief Получение индекса диапазона, которому принадлежит значение аргумента.
    /// @param x Значение аргумента.
    /// @return Индекс диапазона, в котором находится x.
    /// @throws std::logic_error Если значение x не попадает ни в один диапазон.
    size_t get_range_index(double x) const
    {
        // find_if ищет первый элемент, удовлетворяющий условию 
        // (проверил по хелпу https://en.cppreference.com/w/cpp/algorithm/find)
        auto range_coefficients = std::find_if(ranges.begin(), ranges.end(),
            [&](const function_range_t<Coeffs>& r) {
                return x < r.range_end;
            }
        );

        if (range_coefficients == ranges.end())
            throw std::logic_error("approximation range not found");

        if (x > range_coefficients->range_end)
            throw std::logic_error("approximation range not found");

        return range_coefficients - ranges.begin();
    }

    /// @brief Геттер для ranges
    const std::vector<function_range_t<Coeffs>>& get_ranges() const {
        return ranges;
    }

    /// @brief Возвращает полный диапазон - от начала младшего диапазона до конца старшего
    std::pair<double, double> get_whole_range() const {
        if (ranges.empty())
            throw std::runtime_error("cannot get whole range for empty range list");
        
        // пользуемся тем, что диапазоны отсортированы в конструкторе
        return std::make_pair(
            ranges.front().range_start,
            ranges.back().range_end
            );
    }
};

} // namespace fixed_solvers

