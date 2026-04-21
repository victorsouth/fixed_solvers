#define GTEST_BREAK_ON_FAILURE 1
#define GTEST_CATCH_EXCEPTIONS 0
#define GTEST_HAS_SEH 0
#define _VARIADIC_MAX 10 /* for gtest */
#include "gtest/gtest.h"

#include "helpers/cubic_equation_functions.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace {

/// @brief Значение c0 + c1*x + c2*x^2 + c3*x^3.
double eval_poly3(const std::vector<double>& k, double x)
{
    return k[0] + x * (k[1] + x * (k[2] + x * k[3]));
}

/// @brief Число различных вещественных корней по результату solve_cubic_equation
/// (с допуском на совпадение).
int count_unique_real_from_solver(const std::vector<double>& coeffs, double sep_tol)
{
    std::vector<double> roots = fixed_solvers::solve_cubic_equation(coeffs);
    std::sort(roots.begin(), roots.end());
    int n = 0;
    for (size_t i = 0; i < roots.size(); ++i) {
        if (i == 0 || std::fabs(roots[i] - roots[i - 1]) > sep_tol) {
            ++n;
        }
    }
    return n;
}

/// @brief Уникальные корни с слиянием близких значений.
std::vector<double> unique_sorted_roots(const std::vector<double>& roots, double sep_tol)
{
    std::vector<double> u = roots;
    std::sort(u.begin(), u.end());
    std::vector<double> out;
    for (double r : u) {
        if (out.empty() || std::fabs(r - out.back()) > sep_tol) {
            out.push_back(r);
        }
    }
    return out;
}

} // namespace

/// @brief Три различных действительных корня.
TEST(CubicEquationRealRootCount, ReturnsThreeForThreeDistinctRoots)
{
    // Arrange - коэффициенты (x - 1)(x - 2)(x - 3)
    const std::vector<double> p{ -6.0, 11.0, -6.0, 1.0 };
    // Act - подсчёт числа корней по дискриминанту и по списку корней солвера
    const int n_disc = fixed_solvers::count_distinct_real_roots_cubic(p);
    const int n_solver = count_unique_real_from_solver(p, 1e-7);
    // Assert - три различных вещественных корня обоими способами
    ASSERT_EQ(n_disc, 3);
    ASSERT_EQ(n_solver, 3);
}

/// @brief Один действительный корень, два комплексных.
TEST(CubicEquationRealRootCount, ReturnsOneForSingleRealRoot)
{
    // Arrange - полином x^3 - 1
    const std::vector<double> p{ -1.0, 0.0, 0.0, 1.0 };
    // Act - подсчёт корней по дискриминанту и солверу
    const int n_disc = fixed_solvers::count_distinct_real_roots_cubic(p);
    const int n_solver = count_unique_real_from_solver(p, 1e-7);
    // Assert - один вещественный корень
    ASSERT_EQ(n_disc, 1);
    ASSERT_EQ(n_solver, 1);
}

/// @brief Двойной и простой корень: (x - 1)^2 (x + 2) = x^3 - 3*x + 2.
TEST(CubicEquationRealRootCount, ReturnsTwoForDoubleAndSimpleRoot)
{
    // Arrange - (x - 1)^2 (x + 2)
    const std::vector<double> p{ 2.0, -3.0, 0.0, 1.0 };
    // Act - подсчёт корней по дискриминанту и солверу
    const int n_disc = fixed_solvers::count_distinct_real_roots_cubic(p);
    const int n_solver = count_unique_real_from_solver(p, 1e-7);
    // Assert - два различных вещественных корня
    ASSERT_EQ(n_disc, 2);
    ASSERT_EQ(n_solver, 2);
}

/// @brief Тройной корень: (x - 1)^3 = x^3 - 3*x^2 + 3*x - 1.
TEST(CubicEquationRealRootCount, ReturnsOneForTripleRoot)
{
    // Arrange - (x - 1)^3
    const std::vector<double> p{ -1.0, 3.0, -3.0, 1.0 };
    // Act - подсчёт корней по дискриминанту и солверу
    const int n_disc = fixed_solvers::count_distinct_real_roots_cubic(p);
    const int n_solver = count_unique_real_from_solver(p, 1e-7);
    // Assert - один различный вещественный корень (тройной)
    ASSERT_EQ(n_disc, 1);
    ASSERT_EQ(n_solver, 1);
}

/// @brief Масштаб старшего коэффициента не меняет классификацию.
TEST(CubicEquationRealRootCount, IgnoresLeadingCoefficientScale)
{
    // Arrange - 2*(x^3 - 1)
    const std::vector<double> p{ -2.0, 0.0, 0.0, 2.0 };
    // Act - подсчёт числа вещественных корней
    const int n = fixed_solvers::count_distinct_real_roots_cubic(p);
    // Assert - класс как у x^3 - 1: один вещественный корень
    ASSERT_EQ(n, 1);
}

/// @brief Три корня, один из них нуль: x*(x - 1)*(x + 1) = x^3 - x.
TEST(CubicEquationRealRootCount, ReturnsThreeWhenZeroIsRoot)
{
    // Arrange - x (x - 1) (x + 1)
    const std::vector<double> p{ 0.0, -1.0, 0.0, 1.0 };
    // Act - подсчёт корней по дискриминанту и солверу
    const int n_disc = fixed_solvers::count_distinct_real_roots_cubic(p);
    const int n_solver = count_unique_real_from_solver(p, 1e-7);
    // Assert - три различных вещественных корня
    ASSERT_EQ(n_disc, 3);
    ASSERT_EQ(n_solver, 3);
}

/// @brief Неверный размер коэффициентов.
TEST(CubicEquationRealRootCount, ThrowsWhenCoefficientCountIsNotFour)
{
    // Arrange - три коэффициента вместо четырёх
    const std::vector<double> p{ 1.0, 2.0, 3.0 };
    // Act - вызов подсчёта корней с неверной длиной вектора
    // Assert - std::invalid_argument
    ASSERT_THROW(fixed_solvers::count_distinct_real_roots_cubic(p), std::invalid_argument);
}

/// @brief Вырожденная «кубическая» запись со нулевым старшим коэффициентом.
TEST(CubicEquationRealRootCount, ThrowsWhenLeadingCoefficientIsZero)
{
    // Arrange - нулевой старший коэффициент
    const std::vector<double> p{ 1.0, 2.0, 3.0, 0.0 };
    // Act - вызов подсчёта корней
    // Assert - std::invalid_argument
    ASSERT_THROW(fixed_solvers::count_distinct_real_roots_cubic(p), std::invalid_argument);
}

// --- solve_cubic_equation: невязка P(r) ~= 0 ---

/// @brief Три корня: невязка на каждом возвращённом корне.
TEST(SolveCubicEquation, SmallResidualForThreeDistinctRoots)
{
    // Arrange - (x - 1)(x - 2)(x - 3)
    const std::vector<double> p{ -6.0, 11.0, -6.0, 1.0 };
    // Act - поиск вещественных корней
    const std::vector<double> roots = fixed_solvers::solve_cubic_equation(p);
    // Assert - три корня, невязка полинома на каждом
    ASSERT_EQ(roots.size(), 3u);
    for (double r : roots) {
        EXPECT_NEAR(eval_poly3(p, r), 0.0, 1e-9);
    }
}

/// @brief Один вещественный корень: x^3 - 1.
TEST(SolveCubicEquation, SmallResidualForSingleRealRoot)
{
    // Arrange - x^3 - 1
    const std::vector<double> p{ -1.0, 0.0, 0.0, 1.0 };
    // Act - поиск вещественных корней
    const std::vector<double> roots = fixed_solvers::solve_cubic_equation(p);
    // Assert - один корень x = 1 и нулевая невязка
    ASSERT_EQ(roots.size(), 1u);
    EXPECT_NEAR(roots[0], 1.0, 1e-9);
    EXPECT_NEAR(eval_poly3(p, roots[0]), 0.0, 1e-9);
}

/// @brief Пограничный случай S ~= 0: двойной + простой, два значения в ответе.
TEST(SolveCubicEquation, SmallResidualForDoubleRootCase)
{
    // Arrange - (x - 1)^2 (x + 2)
    const std::vector<double> p{ 2.0, -3.0, 0.0, 1.0 };
    // Act - поиск вещественных корней
    const std::vector<double> roots = fixed_solvers::solve_cubic_equation(p);
    // Assert - два значения в ответе, невязка; после слияния близких -2 и 1
    ASSERT_EQ(roots.size(), 2u);
    for (double r : roots) {
        EXPECT_NEAR(eval_poly3(p, r), 0.0, 1e-9);
    }
    auto u = unique_sorted_roots(roots, 1e-7);
    ASSERT_EQ(u.size(), 2u);
    EXPECT_NEAR(u[0], -2.0, 1e-9);
    EXPECT_NEAR(u[1], 1.0, 1e-9);
}

/// @brief Пограничный случай: тройной корень, возможны два совпадающих значения в векторе.
TEST(SolveCubicEquation, SmallResidualForTripleRoot)
{
    // Arrange - (x - 1)^3
    const std::vector<double> p{ -1.0, 3.0, -3.0, 1.0 };
    // Act - поиск вещественных корней
    const std::vector<double> roots = fixed_solvers::solve_cubic_equation(p);
    // Assert - хотя бы один корень, невязка; уникально один x = 1
    ASSERT_GE(roots.size(), 1u);
    for (double r : roots) {
        EXPECT_NEAR(eval_poly3(p, r), 0.0, 1e-9);
    }
    auto u = unique_sorted_roots(roots, 1e-7);
    ASSERT_EQ(u.size(), 1u);
    EXPECT_NEAR(u[0], 1.0, 1e-9);
}

/// @brief Масштаб уравнения не ломает корни.
TEST(SolveCubicEquation, ScaledEquationSameResiduals)
{
    // Arrange - 2*(x^3 - 1)
    const std::vector<double> p{ -2.0, 0.0, 0.0, 2.0 };
    // Act - поиск вещественных корней
    const std::vector<double> roots = fixed_solvers::solve_cubic_equation(p);
    // Assert - один корень x = 1, нулевая невязка
    ASSERT_EQ(roots.size(), 1u);
    EXPECT_NEAR(eval_poly3(p, roots[0]), 0.0, 1e-9);
    EXPECT_NEAR(roots[0], 1.0, 1e-9);
}

/// @brief Три корня включая ноль.
TEST(SolveCubicEquation, ThreeRootsWithZero)
{
    // Arrange - x (x - 1) (x + 1)
    const std::vector<double> p{ 0.0, -1.0, 0.0, 1.0 };
    // Act - поиск вещественных корней
    const std::vector<double> roots = fixed_solvers::solve_cubic_equation(p);
    // Assert - три корня, невязка; уникальные -1, 0, 1
    ASSERT_EQ(roots.size(), 3u);
    for (double r : roots) {
        EXPECT_NEAR(eval_poly3(p, r), 0.0, 1e-9);
    }
    auto u = unique_sorted_roots(roots, 1e-7);
    ASSERT_EQ(u.size(), 3u);
    EXPECT_NEAR(u[0], -1.0, 1e-9);
    EXPECT_NEAR(u[1], 0.0, 1e-9);
    EXPECT_NEAR(u[2], 1.0, 1e-9);
}

/// @brief Согласованность с разложением (x-1)(x-2)(x-3).
TEST(SolveCubicEquation, SortedRootsMatchLinearFactors)
{
    // Arrange - (x - 1)(x - 2)(x - 3)
    const std::vector<double> p{ -6.0, 11.0, -6.0, 1.0 };
    // Act - решение и сортировка корней
    auto roots = fixed_solvers::solve_cubic_equation(p);
    std::sort(roots.begin(), roots.end());
    // Assert - три корня 1, 2, 3
    ASSERT_EQ(roots.size(), 3u);
    EXPECT_NEAR(roots[0], 1.0, 1e-9);
    EXPECT_NEAR(roots[1], 2.0, 1e-9);
    EXPECT_NEAR(roots[2], 3.0, 1e-9);
}
