#pragma once

#define _USE_MATH_DEFINES
#include <algorithm>
#include <cmath>
#include <limits>
#include <math.h>
#include <stdexcept>
#include <vector>


namespace fixed_solvers {

namespace cubic_detail {

/// @brief Абсолютная добавка к порогу сравнения дискриминанта куба с нулём.
inline constexpr double discriminant_zero_tol_absolute = 1e-11;
/// @brief Множитель при scale^3 в пороге для дискриминанта (scale — масштаб коэффициентов).
inline constexpr double discriminant_zero_tol_scale_cubed_coeff = 1e-9;

/// @brief Абсолютная добавка к порогу вырождения производной (a^2 - 3b).
inline constexpr double derivative_degenerate_tol_absolute = 1e-10;
/// @brief Множитель при max(scale^2, 1) в пороге вырождения производной.
inline constexpr double derivative_degenerate_tol_scale_sq_coeff = 1e-8;

/// @brief Коэффициенты монического c + b*x + a*x^2 + x^3 = 0.
struct monic_cubic_coeffs_t {
    /// @brief Коэффициент при x^2.
    double a;
    /// @brief Коэффициент при x.
    double b;
    /// @brief Свободный член.
    double c;
};

/// @brief Пороги для подсчёта числа вещественных корней по дискриминанту.
struct discriminant_count_tols_t {
    /// @brief Порог для сравнения дискриминанта с нулём.
    double tol_disc;
    /// @brief Порог вырождения производной.
    double tol_degen;
};

/// @brief Нормирует коэффициенты c0..c3 к моническому виду (деление на старший).
inline void normalize_to_monic_in_place(std::vector<double>& coeffs)
{
    const double lead = coeffs.back();
    const double lead_rev = 1.0 / lead;
    for (double& coeff : coeffs) {
        coeff *= lead_rev;
    }
}

/// @brief Масштаб max(1, |a|, |b|, |c|) для монического куба.
inline double monic_coefficient_scale(double a, double b, double c)
{
    return std::max({ 1.0, std::fabs(a), std::fabs(b), std::fabs(c) });
}

/// @brief Пороги tol_disc и tol_degen по масштабу монических коэффициентов.
inline discriminant_count_tols_t discriminants_tolerances_for_scale(double scale)
{
    return discriminant_count_tols_t{
        discriminant_zero_tol_absolute
            + discriminant_zero_tol_scale_cubed_coeff * scale * scale * scale,
        derivative_degenerate_tol_absolute
            + derivative_degenerate_tol_scale_sq_coeff * std::max(scale * scale, 1.0),
    };
}

/// @brief Дискриминант монического x^3 + a*x^2 + b*x + c.
inline double monic_discriminant(double a, double b, double c)
{
    return 18.0 * a * b * c - 4.0 * b * b * b + a * a * b * b - 4.0 * a * a * a * c
        - 27.0 * c * c;
}

/// @brief Проверка: ровно 4 коэффициента и ненулевой старший; монические a, b, c.
/// @throws std::invalid_argument Если размер не 4 или старший коэффициент нулевой в смысле eps.
inline monic_cubic_coeffs_t validated_monic_coeffs_from_poly(const std::vector<double>& poly)
{
    if (poly.size() != 4) {
        throw std::invalid_argument("count_distinct_real_roots_cubic: expected 4 coefficients");
    }
    const double lead = poly[3];
    const double lead_abs = std::fabs(lead);
    const double scale_in = std::max(
        { 1.0
        , std::fabs(poly[0])
        , std::fabs(poly[1])
        , std::fabs(poly[2]) });
    if (lead_abs <= std::numeric_limits<double>::epsilon() * scale_in) {
        throw std::invalid_argument("count_distinct_real_roots_cubic: leading coefficient is zero");
    }
    const double lead_rev = 1.0 / lead;
    return monic_cubic_coeffs_t{ poly[2] * lead_rev, poly[1] * lead_rev, poly[0] * lead_rev };
}

/// @brief Число различных действительных корней по моническим коэффициентам.
inline int count_distinct_real_roots_for_monic(const monic_cubic_coeffs_t& m)
{
    const double scale = monic_coefficient_scale(m.a, m.b, m.c);
    const discriminant_count_tols_t tol = discriminants_tolerances_for_scale(scale);
    const double delta = monic_discriminant(m.a, m.b, m.c);

    if (delta > tol.tol_disc) {
        return 3;
    }
    if (delta < -tol.tol_disc) {
        return 1;
    }
    const double pprime_disc = m.a * m.a - 3.0 * m.b;
    if (std::fabs(pprime_disc) <= tol.tol_degen) {
        return 1;
    }
    return 2;
}

/// @brief Ветвь S ~= 0: кратные корни (два значения в возвращаемом векторе).
inline std::vector<double> solve_monic_s_near_zero(double a, double R)
{
    using std::cbrt;
    return { -2. * cbrt(R) - a / 3., cbrt(R) - a / 3. };
}

/// @brief Ветвь S > 0: три различных действительных корня (формула Виета).
inline std::vector<double> solve_monic_three_real_roots(double a, double Q, double R)
{
    using std::acos;
    using std::cos;
    using std::pow;
    using std::sqrt;

    const double phi = 1.0 / 3.0 * acos(R / pow(Q, 1.5));
    const double shift = -a / 3;
    const double radius = -2 * sqrt(Q);
    return {
        radius * cos(phi) + shift,
        radius * cos(phi + 2.0 / 3.0 * M_PI) + shift,
        radius * cos(phi - 2.0 / 3.0 * M_PI) + shift,
    };
}

/// @brief Скаляр Кардано (p/3)^3 + (q/2)^2 для монического x^3 + a*x^2 + b*x + c = 0
/// после перехода к приведённому виду (те же p_dep, q_dep, что в ветке одного вещественного корня).
inline double monic_depressed_cubic_cardano_Q(double a, double b, double c)
{
    const double p_dep = b - a * a / 3.0;
    const double q_dep = c - a * b / 3.0 + 2.0 * a * a * a / 27.0;
    return std::pow(p_dep / 3.0, 3) + std::pow(q_dep / 2.0, 2);
}

/// @brief Ветвь S < 0: один действительный корень (депрессивный куб + Кардано).
inline std::vector<double> solve_monic_one_real_depressed_cardano(double a, double b, double c)
{
    const double p_dep = b - a * a / 3.0;
    const double q_dep = c - a * b / 3.0 + 2.0 * a * a * a / 27.0;
    const double Q_card = std::pow(p_dep / 3.0, 3) + std::pow(q_dep / 2.0, 2);
    const double rad = std::sqrt(std::max(0.0, Q_card));
    const double d1 = -q_dep / 2.0 + rad;
    const double d2 = -q_dep / 2.0 - rad;
    const double y = std::cbrt(d1) + std::cbrt(d2);
    return { y - a / 3.0 };
}

/// @brief Действительные корни монического c + b*x + a*x^2 + x^3 = 0.
inline std::vector<double> solve_monic(double a, double b, double c)
{
    using std::fabs;
    using std::pow;

    const double Q = (a * a - 3 * b) / 9;
    const double R = (2 * pow(a, 3) - 9 * a * b + 27 * c) / 54;
    const double S = pow(Q, 3) - R * R;

    if (fabs(S) < std::numeric_limits<double>::epsilon()) {
        return solve_monic_s_near_zero(a, R);
    }
    if (S > 0) {
        return solve_monic_three_real_roots(a, Q, R);
    }
    return solve_monic_one_real_depressed_cardano(a, b, c);
}

/// @brief Число различных действительных корней (полином степени 3, см. публичный API).
inline int count_distinct_real_roots(const std::vector<double>& poly_coeffs)
{
    const monic_cubic_coeffs_t m = validated_monic_coeffs_from_poly(poly_coeffs);
    return count_distinct_real_roots_for_monic(m);
}

} // namespace cubic_detail

/// @brief Число различных действительных корней кубического уравнения.
/// @details Полином c0 + c1*x + c2*x^2 + c3*x^3 = 0, индекс коэффициента равен степени
/// (как у solve_cubic_equation). После приведения к моническому виду x^3 + a*x^2 + b*x + c
/// используется дискриминант; при нулевом дискриминанте тройной корень отличается от
/// «двойной + простой» по условию a^2 = 3b на производной.
/// @param poly_coeffs Ровно четыре коэффициента, полином третьей степени (c3 != 0).
/// @return 1, 2 или 3 различных действительных корня.
/// @throws std::invalid_argument Если размер вектора не 4 или старший коэффициент нулевой.
inline int count_distinct_real_roots_cubic(const std::vector<double>& poly_coeffs)
{
    return cubic_detail::count_distinct_real_roots(poly_coeffs);
}

/// @brief Действительные корни кубического уравнения (комплексные не возвращаются).
/// @details Три вещественных корня — тригонометрическая формула Виета; один вещественный —
/// явный Кардано по депрессивному виду (как Maple roots1_eqZ / roots1_eqZ_mix).
/// @param poly_coeffs Четыре коэффициента c0..c3 при c0 + c1*x + c2*x^2 + c3*x^3 = 0 (индекс = степень).
/// @return Все действительные корни на вещественной оси (1, 2 или 3 значения).
inline std::vector<double> solve_cubic_equation(const std::vector<double>& poly_coeffs)
{
    std::vector<double> coeffs = poly_coeffs;
    cubic_detail::normalize_to_monic_in_place(coeffs);
    const double a = coeffs[2];
    const double b = coeffs[1];
    const double c = coeffs[0];
    return cubic_detail::solve_monic(a, b, c);
}

/// @brief Скаляр Кардано (p/3)^3 + (q/2)^2 для монического x^3 + a*x^2 + b*x + c = 0.
inline double monic_cubic_cardano_Q(double a, double b, double c)
{
    return cubic_detail::monic_depressed_cubic_cardano_Q(a, b, c);
}

} // namespace fixed_solvers
