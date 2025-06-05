/*!
* \file
* \brief В данном .h файле реализовано описание систем различной
размерности для решателя Ньютона - Рафсона
*
* \author Автор файла - В.В. Южанин, Автор документации - И.Б. Цехместрук
* \date Дата докуменатции - 2023-04-06
*
* Документация к этому файлу находится в:
* 1. 01. Антуан Рауль-Дальтон\02. Документы - черновики\Иван\01. Описание численного метода
*/

#pragma once

#include <numeric>

#include <Eigen/Dense>
#include <Eigen/Sparse>


using std::accumulate;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
typedef Eigen::Triplet<double> triplet_t;

/// @brief Расчет приращения для численного расчета производной на основе относительного отклонения
/// @param value Точка, где вычисляется производная
/// @param epsilon Относительное отклонение
/// @return Приращение
inline double numeric_derivative_delta(double value, double epsilon)
{
    using std::max;
    return epsilon * max(1.0, std::abs(value));
}


/// @brief Численный расчет производной от функции одного аргумента по двусторонней формуле
/// @tparam Function Тип функции
/// @param f Функция
/// @param value Точка, в которой вычисляется произвожная
/// @param epsilon Относительное приращение (см. numeric_derivative_delta)
/// @return Производная 
template <typename Function>
inline auto two_sided_derivative(Function f, double value, double epsilon) {
    double dx = numeric_derivative_delta(value, epsilon);

    typedef std::invoke_result_t<Function, double> ResultType;
    ResultType f_plus = f(value + dx);
    ResultType f_minus = f(value - dx);
    ResultType difference = f_plus - f_minus;
    ResultType result = difference / (2 * dx);
    return result;
}

/// @brief Типы данных для системы уравнений
template <std::ptrdiff_t Dimension>
struct fixed_system_types;

/// @brief Специализация вырожденного случая системы уравнений - одного уравнения, скалярный случай
template <>
struct fixed_system_types<1> {
    /// @brief Тип переменной
    typedef double var_type;
    /// @brief Тип правой части
    typedef double right_party_type;
    /// @brief Тип матрицы (линейная система, якобиан)
    typedef double matrix_type;
    /// @brief Тип коэффициентов уравнения
    typedef double equation_coeffs_type;

    /// @brief Инициализация неизвестной переменной по умолчанию для скалярного случая
    static var_type default_var(double value = std::numeric_limits<double>::quiet_NaN())
    {
        return value;
    }
};

/// @brief Специализация случая переменной размерности
/// @details В данный момент не реализовано и не используется
template <>
struct fixed_system_types<-1> {
    /// @brief Тип переменной
    typedef VectorXd var_type;
    /// @brief Тип правой части
    typedef VectorXd right_party_type;
    /// @brief Тип матрицы (линейная система, якобиан)
    typedef MatrixXd matrix_type;
    /// @brief Тип коэффициентов уравнения
    typedef MatrixXd equation_coeffs_type;
    /// @brief Тип разреженной матрицы
    typedef std::vector<Eigen::Triplet<double>> sparse_matrix_type;
    /// @brief Инициализация неизвестной переменной по умолчанию для скалярного случая
    static var_type default_var(
        double value = std::numeric_limits<double>::quiet_NaN(),
        size_t dimension = 0)
    {
        VectorXd result = VectorXd::Ones(dimension);
        result *= value;
        return result;
    }
};

/// @brief Общий случай системы уравнений фиксированной размерности
template <std::ptrdiff_t Dimension>
struct fixed_system_types {
    /// @brief Тип переменной
    typedef array<double, Dimension> var_type;
    /// @brief Тип правой части
    typedef var_type right_party_type;
    /// @brief Тип матрицы (линейная система, якобиан)
    typedef array<var_type, Dimension> matrix_type;
    /// @brief Тип коэффициентов уравнения
    typedef matrix_type equation_coeffs_type;
    /// @brief Инициализация неизвестной переменной по умолчанию для скалярного случая
    static var_type default_var(double value = std::numeric_limits<double>::quiet_NaN())
    {
        auto getter = [&](size_t index) { return value; };
        return create_array<Dimension>(getter);
    }
};


template <std::ptrdiff_t Dimension>
class fixed_system_t;

/// @brief Система уравнений с переменной размерности
template <>
class fixed_system_t<-1> {
protected:
    /// @brief Относительное (!) приращение для расчета производных
    double epsilon{ 1e-6 };
public:
    /// @brief Тип разреженной матрицы
    typedef typename fixed_system_types<-1>::sparse_matrix_type sparse_matrix_type;
public:
    /// @brief Расчет целевой функции по невязкам
    virtual double objective_function(const VectorXd& r) const
    {
        return r.squaredNorm();
    }
    /// @brief Расчет целевой функции по аргументу
    double operator()(const VectorXd& x) {
        auto r = residuals(x);
        return objective_function(r);
    }
    /// @brief Невязки системы уравнений
    virtual VectorXd residuals(const VectorXd& x) = 0;
    /// @brief Якобиан системы уравнений
    virtual sparse_matrix_type jacobian_sparse(const VectorXd& x) {
        return jacobian_sparse_numeric(x);
    }
    /// @brief Колонка разреженного якобиана 
    /// Очень примитивная реализация через исходный разреженный якобиан
    virtual sparse_matrix_type jacobian_sparse_column(const VectorXd& reduced, size_t desired_col_index) {
        std::vector<Eigen::Triplet<double>> J = jacobian_sparse(reduced);
        std::vector<Eigen::Triplet<double>> result;
        for (const Eigen::Triplet<double>& triplet : J) {
            Eigen::Triplet<double> new_triplet(triplet.row(), 0, triplet.value());
            size_t col_index = triplet.col();
            if (desired_col_index == col_index) {
                result.push_back(new_triplet);
            }
        }

        return result;
    }
    /// @brief Специфический критерий успешного завершения расчета
    /// @param r Текущее значения невязок
    /// @param x Текущее значение аргумента
    /// @return Флаг успешного завершения
    virtual bool custom_success_criteria(const VectorXd& r, const VectorXd& x)
    {
        return false;
    }
protected:
    /// @brief Численный расчет разреженного якобиана
    sparse_matrix_type jacobian_sparse_numeric(const VectorXd& x) {
        std::vector<Eigen::Triplet<double>> result;
        VectorXd arg = x;

        for (int component = 0; component < x.size(); ++component) {
            double e = numeric_derivative_delta(arg[component], epsilon);
            arg[component] = x[component] + e;
            auto f_plus = residuals(arg);
            arg[component] = x[component] - e;
            auto f_minus = residuals(arg);
            arg[component] = x[component];

            VectorXd Jcol = (f_plus - f_minus) / (2 * e);
            for (int row = 0; row < x.size(); ++row) {
                result.emplace_back(row, component, Jcol(row));
            }
        }
        return result;
    }
};


/// @brief Система алгебраических уравнений - базовый класс
template <std::ptrdiff_t Dimension>
class fixed_system_t {
protected:
    /// @brief Относительное (!) приращение для расчета производных
    double epsilon{ 1e-6 };
public:
    /// @brief Тип переменной
    typedef typename fixed_system_types<Dimension>::var_type var_type;
    /// @brief Тип функции
    typedef typename fixed_system_types<Dimension>::right_party_type function_type;
    /// @brief Тип матрицы (якобиан)
    typedef typename fixed_system_types<Dimension>::equation_coeffs_type matrix_value;
    /// @brief Тип разреженной матрицы (якобиан)
    typedef typename fixed_system_types<-1>::sparse_matrix_type sparse_matrix_type;
public:
    /// @brief Расчет целевой функции по невязкам
    virtual double objective_function(const var_type& r) const;
    /// @brief Расчет целевой функции по аргументу
    double operator()(const var_type& x) {
        auto r = residuals(x);
        return objective_function(r);
    }
    /// @brief Невязки системы уравнений
    virtual function_type residuals(const var_type& x) = 0;
    /// @brief Якобиан системы уравнений
    virtual matrix_value jacobian_dense(const var_type& x) {
        return jacobian_dense_numeric(x);
    }
    /// @brief Якобиан системы уравнений
    virtual sparse_matrix_type jacobian_sparse(const var_type& x) {
        return jacobian_sparse_numeric(x);
    }
    /// @brief Колонка разреженного якобиана 
    virtual function_type jacobian_column(const var_type& x, size_t desired_col_index) {
        if constexpr (Dimension == 1) {
            throw std::runtime_error("Jacobian column for dimension = 1 is sensless");
        }
        else {
            var_type arg = x;

            double e = numeric_derivative_delta(arg[desired_col_index], epsilon);
            arg[desired_col_index] = x[desired_col_index] + e;
            function_type f_plus = residuals(arg);
            arg[desired_col_index] = x[desired_col_index] - e;
            function_type f_minus = residuals(arg);
            arg[desired_col_index] = x[desired_col_index];

            function_type Jcol = (f_plus - f_minus) / (2 * e);
            return Jcol;
        }
    }

    /// @brief Специфический критерий успешного завершения расчета
    /// @param r Текущее значения невязок
    /// @param x Текущее значение аргумента
    /// @return Флаг успешного завершения
    virtual bool custom_success_criteria(const var_type& r, const var_type& x)
    {
        return false;
    }

protected:
    /// @brief Численный расчет плотного якобиана
    matrix_value jacobian_dense_numeric(const var_type& x) {
        if constexpr (Dimension == 1) {
            auto f = [&](double x) { return residuals(x); };
            double result = two_sided_derivative(f, x, epsilon);
            return result;
        }
        else {
            var_type arg = x;

            matrix_value J;

            for (int component = 0; component < x.size(); ++component) {
                double e = numeric_derivative_delta(arg[component], epsilon);
                arg[component] = x[component] + e;
                function_type f_plus = residuals(arg);
                arg[component] = x[component] - e;
                function_type f_minus = residuals(arg);
                arg[component] = x[component];

                function_type Jcol = (f_plus - f_minus) / (2 * e);
                for (size_t row = 0; row < static_cast<size_t>(x.size()); ++row) {
                    J[row][component] = Jcol[row];
                }
            }
            return J;

        }
    }
    /// @brief Численный расчет разреженного якобиана
    sparse_matrix_type jacobian_sparse_numeric(const var_type& x) {

        if constexpr (Dimension > 1) {
            std::vector<Eigen::Triplet<double>> result;
            var_type arg = x;

            for (int component = 0; component < Dimension; ++component) {
                double e = numeric_derivative_delta(arg[component], epsilon);
                arg[component] = x[component] + e;
                auto f_plus = residuals(arg);
                arg[component] = x[component] - e;
                auto f_minus = residuals(arg);
                arg[component] = x[component];

                var_type Jcol = (f_plus - f_minus) / (2 * e);
                for (int row = 0; row < static_cast<int>(Dimension); ++row) {
                    result.emplace_back(row, component, Jcol[row]);
                }
            }
            return result;

        }
        else {
            throw std::runtime_error("Must not be called");
        }

    }
};

/// @brief Расчет целевой функции для скалярного случая (сумма квадратов)
template <>
inline double fixed_system_t<1>::objective_function(const double& r) const
{
    return r * r;
}


/// @brief Расчет целевой функции для векторного случая (сумма квадратов)
template <std::ptrdiff_t Dimension>
inline double fixed_system_t<Dimension>::objective_function(const var_type& r) const
{
    double result = std::accumulate(
        r.begin(), r.end(), 0.0, 
        [](double accum, double value) { return accum + value * value; });

    return result;
}

/// @brief Оборачивает лямбда-функцию в fixed_system_t<1>
/// Реализация якобина - численная из fixed_system_t<1>
/// @tparam Lambda Тип лямба-функции
template <typename Lambda>
class fixed_scalar_wrapper_t : public fixed_system_t<1>
{
protected:
    /// @brief Указатель на лямбду
    Lambda function;
public:
    /// @brief Конструктор позволяет переопределить 
    /// @param function 
    /// @param epsilon 
    fixed_scalar_wrapper_t(Lambda function, 
        double epsilon = std::numeric_limits < double >::quiet_NaN())
        : function(function)
    {
        if (!std::isnan(epsilon)) {
            this->epsilon = epsilon;
        }
    }
    /// @brief Невязки уравнения
    virtual function_type residuals(const var_type& x) override
    {
        return function(x);
    }

};
