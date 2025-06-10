#pragma once

using std::isfinite;

/// @brief Решение линейного уравнения
/// ax = b
inline double solve_linear_system(
    const double& a, const double& b)
{
    double result = b / a;
    if (!std::isfinite(result)) {
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
inline std::array<double, 2> solve_linear_system(
    const std::array<std::array<double, 2>, 2>& A, const std::array<double, 2>& b)
{

    double d = A[0][0] * A[1][1] - A[1][0] * A[0][1];
    double d1 = b[0] * A[1][1] - b[1] * A[0][1];
    double d2 = A[0][0] * b[1] - A[1][0] * b[0];

    std::array<double, 2> result{ d1 / d, d2 / d };

    if (!std::isfinite(result[0]) || !std::isfinite(result[1])) {
        throw std::logic_error("infinite value");
    }

    return result;
}


inline double determinant(
    double a11, double a12, double a13, 
    double a21, double a22, double a23,
    double a31, double a32, double a33
    )
{
    return 
        + a11* a22* a33 
        - a11 * a23 * a32 
        - a12 * a21 * a33 
        + a12 * a23 * a31 
        + a13 * a21 * a32 
        - a13 * a22 * a31;
}

/// @brief Решение СЛАУ 3х3 методом Крамера
/// Ax = b
/// https://ru.wikipedia.org/wiki/%D0%9C%D0%B5%D1%82%D0%BE%D0%B4_%D0%9A%D1%80%D0%B0%D0%BC%D0%B5%D1%80%D0%B0
inline std::array<double, 3> solve_linear_system(
    const std::array<std::array<double, 3>, 3>& A, const std::array<double, 3>& b)
{
    double d = determinant(
        A[0][0], A[0][1], A[0][2],
        A[1][0], A[1][1], A[1][2],
        A[2][0], A[2][1], A[2][2]);

    double d1 = determinant(
        b[0], A[0][1], A[0][2],
        b[1], A[1][1], A[1][2],
        b[2], A[2][1], A[2][2]);

    double d2 = determinant(
        A[0][0], b[0], A[0][2],
        A[1][0], b[1], A[1][2],
        A[2][0], b[2], A[2][2]);

    double d3 = determinant(
        A[0][0], A[0][1], b[0],
        A[1][0], A[1][1], b[1],
        A[2][0], A[2][1], b[2]);


    std::array<double, 3> result{ d1 / d, d2 / d, d3 / d };

    if (!std::isfinite(result[0]) || !std::isfinite(result[1]) || !std::isfinite(result[2])) {
        throw std::logic_error("infinite value");
    }

    return result;
}
