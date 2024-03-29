﻿/*!
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

/// @brief Типы данных для системы уравнений
template <std::ptrdiff_t Dimension>
struct fixed_system_types;

/// @brief Специализация вырожденного случая системы уравнений - одного уравнения, скалярный случай
template <>
struct fixed_system_types<1> {
    typedef double var_type;
    typedef double right_party_type;
    typedef double matrix_type;
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
    typedef VectorXd var_type;
    typedef VectorXd right_party_type;
    typedef MatrixXd matrix_type;
    typedef MatrixXd equation_coeffs_type;

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
    typedef array<double, Dimension> var_type;
    typedef var_type right_party_type;
    typedef array<var_type, Dimension> matrix_type;
    typedef matrix_type equation_coeffs_type;

    /// @brief Инициализация неизвестной переменной по умолчанию для скалярного случая
    static var_type default_var(double value = std::numeric_limits<double>::quiet_NaN())
    {
        auto getter = [&](size_t index) { return value; };
        return create_array<Dimension>(getter);
    }

};


/// @brief Система алгебраических уравнений - базовый класс
template <std::ptrdiff_t Dimension>
class fixed_system_t {
protected:
    /// @brief Относительное (!) приращение для расчета производных
    double epsilon{ 1e-6 };
public:
    typedef typename fixed_system_types<Dimension>::var_type var_type;
    typedef typename fixed_system_types<Dimension>::right_party_type function_type;
    typedef typename fixed_system_types<Dimension>::equation_coeffs_type matrix_value;
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
    /// @brief Специфический критерий успешного завершения расчета
    /// @param r Текущее значения невязок
    /// @param x Текущее значение аргумента
    /// @return Флаг успешного завершения
    virtual bool custom_success_criteria(const var_type& r, const var_type& x)
    {
        return true;
    }

protected:
    matrix_value jacobian_dense_numeric(const var_type& x);

};

/// @brief Численный расчет Якобиана методом двусторонней разности для скалярного случая
/// @param x Текущее значение аргумента
/// @return Значение Якобиана
template <>
inline double fixed_system_t<1>::jacobian_dense_numeric(const double& x)
{
    double e = epsilon * std::max(1.0, abs(x));
    function_type f_plus = residuals(x + e);
    function_type f_minus = residuals(x - e);
    function_type J = (f_plus - f_minus) / (2 * e);
    return J;

}


template <>
inline fixed_system_t<-1>::matrix_value 
    fixed_system_t<-1>::jacobian_dense_numeric(const fixed_system_t<-1>::var_type& x)
{
    throw std::logic_error("var system must have sparse jacobian");

}

/// @brief Численный расчет Якобиана методом двусторонней разности для векторного случая
/// @tparam Dimension Размерность решаемой системы уравнений
/// @param x Текущее значение аргумента
/// @return Значение Якобиана
template <std::ptrdiff_t Dimension>
inline typename fixed_system_t<Dimension>::matrix_value fixed_system_t<Dimension>::jacobian_dense_numeric(const var_type& x)
{
    var_type arg = x;

    matrix_value J;

    for (int component = 0; component < x.size(); ++component) {
        double e = epsilon * std::max(1.0, abs(arg[component]));
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

/// @brief Расчет целевой функции для скалярного случая (сумма квадратов)
template <>
inline double fixed_system_t<1>::objective_function(const double& r) const
{
    return r * r;
}


template <>
inline double fixed_system_t<-1>::objective_function(const VectorXd& r) const
{
    return r.squaredNorm();
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

/// @brief Реализация residuals и jacobian для лямбда-функций
/// @tparam Lambda Лямба-функция
template <typename Lambda>
class fixed_scalar_wrapper_t : public fixed_system_t<1>
{
protected:
    Lambda function;
public:
    fixed_scalar_wrapper_t(
        Lambda function, double epsilon = std::numeric_limits < double >::quiet_NaN())
        :function(function)
    {
        if (!std::isnan(epsilon)) {
            this->epsilon = epsilon;
        }
    }
    virtual function_type residuals(const var_type& x) override
    {
        return function(x);
    }

};
