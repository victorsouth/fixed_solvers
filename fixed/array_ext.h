#pragma once

#include <array>
using std::array;


template<typename T, std::size_t N, int Idx = N>
struct array_maker {
    template<typename... Ts>
    static std::array<T, N> make_array(const T& v, Ts...tail) {
        return array_maker<T, N, Idx - 1>::make_array(v, v, tail...);
    }
};

template<typename T, std::size_t N>
struct array_maker<T, N, 1> {
    template<typename... Ts>
    static std::array<T, N> make_array(const T& v, Ts... tail) {
        return std::array<T, N>{v, tail...};
    }
};

template<typename T>
struct array_maker<T, 0, -1> {
    template<typename... Ts>
    static std::array<T, 0> make_array(const T& v, Ts... tail) {
        return std::array<T, 0>{};
    }
};



template <typename DataType, size_t Dimension>
inline array<DataType, Dimension> operator - (
    const array<DataType, Dimension>& v1, const array<DataType, Dimension>& v2)
{
    array<DataType, Dimension> result(v1);
    for (size_t index = 0; index < v1.size(); ++index) {
        result[index] -= v2[index];
    }
    return result;
}

template <typename DataType, size_t Dimension>
inline array<DataType, Dimension> operator + (
    const array<DataType, Dimension>& v1, const array<DataType, Dimension>& v2)
{
    array<DataType, Dimension> result(v1);
    for (size_t index = 0; index < v1.size(); ++index) {
        result[index] += v2[index];
    }
    return result;
}

template <typename DataType, size_t Dimension>
inline array<DataType, Dimension> operator * (
    double scalar, const array<DataType, Dimension>& v1)
{
    array<DataType, Dimension> result(v1);
    for (size_t index = 0; index < v1.size(); ++index) {
        result[index] *= scalar;
    }
    return result;
}

/// @brief Скалярное произведение векторов array<N>
template <typename DataType, size_t Dimension>
inline constexpr DataType inner_prod(
    const array<DataType, Dimension>& v1,
    const array<DataType, Dimension>& v2)
{
    DataType result = v1[0] * v2[0];
    for (size_t index = 1; index < Dimension; ++index) {
        result += v1[index] * v2[index];
    }
    return result;
}

/// @brief Перегрузка, когда вектор обращается в скаляр
inline constexpr double inner_prod(
    double v1,
    double v2)
{
    return v1 * v2;
}

/// @brief Оператор вычитания векторов фиксированной размерности
template <typename DataType, size_t Dimension>
array<DataType, Dimension> operator - (
    array<DataType, Dimension> v)
{
    for (size_t index = 0; index < v.size(); ++index) {
        v[index] = -v[index];
    }
    return v;
}

/// @brief Оператор вычитания векторов фиксированной размерности
template <typename DataType, size_t Dimension>
array<array<DataType, Dimension>, Dimension> operator - (
    array<array<DataType, Dimension>, Dimension> m)
{
    for (size_t index = 0; index < m.size(); ++index) {
        m[index] = -m[index];
    }
    return m;
}

/// @brief Оператор умножения векторов фиксированной размерности
template <typename DataType, size_t Dimension>
array<DataType, Dimension> operator * (
    const array<DataType, Dimension>& v1, double scalar)
{
    return scalar * v1;
}

/// @brief Оператор поэлементного деления векторов фиксированной размерности
template <typename DataType, size_t Dimension>
array<DataType, Dimension> operator / (
    const array<DataType, Dimension>& v1, double scalar)
{
    auto result = v1;
    for (auto& value : result) {
        value /= scalar;
    }
    return result;
}



/// @brief Умножение матрицы на вектор-столбец справа
/// @param m Матрица Row-major
/// @param v Вектор-столбец
/// @return Вектор-столбец
template <typename DataType, size_t Dimension>
array<DataType, Dimension> operator * (
    const array<array<DataType, Dimension>, Dimension>& m, const array<DataType, Dimension>& v)
{
    array<DataType, Dimension> result;
    for (size_t index = 0; index < m.size(); ++index) {
        result[index] = inner_prod(m[index], v);
    }
    return result;
}

/// @brief Оператор сокращенного сложения для векторов фиксированной размерности
template <typename DataType, size_t Dimension>
array<DataType, Dimension>& operator += (array<DataType, Dimension>& v1, const array<DataType, Dimension>& v2) {
    for (size_t index = 0; index < v1.size(); ++index) {
        v1[index] += v2[index];
    }
    return v1;
}

/// @brief Для create_array
template<int I, size_t N, typename DataType, typename Getter, typename... Tp>
static inline typename std::enable_if<I == -1, array<DataType, N>>::type
create_array_recursion(const Getter& getter, Tp...tail)
{
    array<DataType, N> r{ tail... };
    return r;
}

/// @brief Для create_array
template<int I, size_t N, typename DataType, typename Getter, typename... Tp>
static inline typename std::enable_if < I >= 0, array<DataType, N>>::type
create_array_recursion(const Getter& getter, Tp...tail)
{
    return create_array_recursion<I - 1, N, DataType>(getter, getter(I), tail...);
}

/// @brief Создает массив заданного размера по геттеру
/// Геттер вызывается с параметров индекса элемента
template <size_t Dimension, typename Getter>
static inline auto create_array(const Getter& getter)
{
    typedef typename std::invoke_result_t<Getter, int> DataType;
    return create_array_recursion<Dimension - 1, Dimension, DataType>(getter);
}

/// @brief Реализует поведение ссылочного массива
template <typename T, size_t Dimension>
class array_ref {
private:
    /// @brief Массив указателей
    array<T*, Dimension> refs;
public:
    /// @brief Очевидный конструктор
    array_ref(array<T*, Dimension>& refs)
        : refs(refs)
    {

    }
    /// @brief Создаем ссылочный массив для передаваемого массива
    array_ref(array<T, Dimension>& a)
        : refs(create_array<Dimension>([&](int dimension) { return &a[dimension]; }))
    {

    }
    /// @brief Что собой представляет getter? 
    template <typename Getter>
    array_ref(const Getter& getter)
        : refs(create_array<Dimension>(getter))
    {

    }
    /// @brief Где используется?
    operator array<T, Dimension>() {
        array<T, Dimension> result;
        for (size_t index = 0; index < Dimension; ++index) {
            result[index] = *refs[index];
        }
        return result;
    }

    /// @brief Обратиться к элементу ссылочного массива
    T& operator[](size_t index) {
        return *refs[index];
    }
    
    /// @brief Переприсвоить элементы ссылочного массива
    void operator= (const array<T, Dimension>& other)
    {
        for (size_t index = 0; index < Dimension; ++index) {
            *refs[index] = other[index];
        }
    }
    /// @brief Переприсвоить элементы ссылочного массива (скалярный случай)
    void operator= (double other)
    {
        for (size_t index = 0; index < Dimension; ++index) {
            *refs[index] = other;
        }
    }
    /// @brief Где используется?
    operator T () const {
        return *(refs[0]);
    }
};


