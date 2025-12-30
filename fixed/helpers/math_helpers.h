#pragma once



namespace fixed_solvers {
;

template <typename T> inline int pseudo_sgn(T val) {
    return 2 * (val >= T(0)) - 1;
}

inline Eigen::VectorXd invVectorXd(const Eigen::VectorXd& _v)
{
    Eigen::VectorXd v = _v;
    for (auto index = 0; index < v.size(); ++index) {
        v(index) = 1.0 / v(index);
    }

    return v;
}

/// @brief Возвращает знак числа
/// @return для отрицательных -1, для положительных +1
template <typename T> inline int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

/// @brief Возведение вещественной переменной в квадрат
inline double sqr(double x) {
    return x * x;
}

/// @brief Вычисление квадратного корня с учетом знака подкоренного выражения
/// @param x - подкоренное выражение
/// @return квадратный корень (для отрицательных возвращает -sqrt(-x))
inline double ssqrt(double x) {
    if (x >= 0)
        return sqrt(x);
    else
        return -sqrt(-x);
}

/// @brief Возведение в квадрат с учетом знака
/// @param x - значение
/// @return квадрат значения (для отрицательных возвращает -sqr(x))
inline double ssqr(double x) {
    if (x >= 0)
        return sqr(x);
    else
        return -sqr(x);
}

template <typename Coefficients>
double polyval(const Coefficients& poly_coeffs, double x)
{
    // НЕ используется формат Matlab!
    // индекс коэффициента является показателем степени

    const size_t n = poly_coeffs.size();

    if (n == 0)
        return 0;

    double result = 0;
    double power = 1;
    for (size_t index = 0; index < n; ++index) {
        result += power * poly_coeffs[index];
        power *= x;
    }
    return result;
}

inline bool has_not_finite(const double value)
{
    return !std::isfinite(value);
}

inline bool has_not_finite(const Eigen::VectorXd& v)
{
    for (int index = 0; index < v.size(); ++index) {
        if (!std::isfinite(v(index))) {
            return true;
        }
    }
    return false;
}

}

