#pragma once

/// @brief Решение линейного уравнения
/// ax = b
inline double solve_linear_system(
    const double& a, const double& b)
{
    double result = b / a;
    if (!isfinite(result)) {
        throw std::logic_error("infinite value");
    }
    return result;
}

/// @brief Решение линейного уравнения
/// ax = b
/// @param equation a = first, b = second
inline double solve_linear_system(
    const std::pair<double, double>& equation)
{
    return solve_linear_system(equation.first, equation.second);
}



/// @brief Решение СЛАУ 2х2 методом Крамера
/// Ax = b
/// https://ru.wikipedia.org/wiki/%D0%9C%D0%B5%D1%82%D0%BE%D0%B4_%D0%9A%D1%80%D0%B0%D0%BC%D0%B5%D1%80%D0%B0
inline array<double, 2> solve_linear_system(
    const array<array<double, 2>, 2>& A, const array<double, 2>& b)
{

    double d = A[0][0] * A[1][1] - A[1][0] * A[0][1];
    double d1 = b[0] * A[1][1] - b[1] * A[0][1];
    double d2 = A[0][0] * b[1] - A[1][0] * b[0];

    array<double, 2> result{ d1 / d, d2 / d };

    if (!isfinite(result[0]) || !isfinite(result[1])) {
        throw std::logic_error("infinite value");
    }

    return result;
}
