#pragma once

#ifndef __FIXED_H__
#define __FIXED_H__

#include <numeric>
#include <Eigen/Dense>
#include <Eigen/Sparse>


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

// Действительные корни уравнения 3-й степени



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


inline std::string wide2string(const std::wstring& str)
{
    std::string ascii_string(str.size(), '\0');
    std::transform(str.begin(), str.end(), ascii_string.begin(),
        [](wchar_t s) { return static_cast<char>(s); });
    return ascii_string;
}
inline std::wstring string2wide(const std::string& str)
{
#ifndef __GNUC__
    try {
        std::wstring wstr(str.size(), 0);
        std::use_facet<std::ctype<wchar_t> >(std::locale("rus_rus.1251")).widen
        (&str[0], &str[0] + str.size(), &wstr[0]);
        return wstr;
    }
    catch (...)
    {
        std::wstring wide_string(str.begin(), str.end());
        return wide_string;
    }
#else
    std::wstring wstr = UTF8_to_wchar(str.c_str());
    return wstr;
#endif
}


template<class T>
inline std::string int2str(T i)
{
    std::stringstream ss;
    ss << i;
    return ss.str();
}

template<class T>
inline std::wstring int2wstr(T i)
{
    std::wstringstream ss;
    ss << i;
    return ss.str();
}

inline void string_replace(std::string& str, const std::string& from, const std::string& to)
{
    if (from.empty())
        return;
    size_t start_pos = 0;
    while ((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
    }
}


}

#include <iomanip>
#include "fixed/qp/qp_wrapper.h"
#include "fixed/array_ext.h"
#include "fixed/fixed_system.h"
#include "fixed/fixed_linear_solver.h"
#include "fixed/fixed_constraints.h"
#include "fixed/fixed_nonlinear_solver.h"
#include "fixed/fixed_optimizer.h"

#endif
