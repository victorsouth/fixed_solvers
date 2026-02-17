#pragma once



namespace fixed_solvers {
;

/// @brief Булевый тип с возможностью не задавать значение
enum class nullable_bool_t {
    True = 0b1, False = 0b0, Undefined = 0b10
};

/// @brief Преобразование bool в nullable bool
inline nullable_bool_t to_nullable_bool(bool value)
{
    return value
        ? nullable_bool_t::True
        : nullable_bool_t::False;
}

/// @brief Вывод в поток для nullable bool
inline std::ostream& operator<<(std::ostream& o, nullable_bool_t x) {
	switch (x)
	{
	case nullable_bool_t::True:
		o << "True";
		break;
	case nullable_bool_t::False:
		o << "False";
		break;
	case nullable_bool_t::Undefined:
		o << "Undefined";
		break;
	default:
		o << "unknown nullable bool value";
		break;
	}
	return o;
}


/// @brief Вывод в widestring-поток для nullable bool
inline std::wostream& operator<<(std::wostream& o, nullable_bool_t x) {
	switch (x)
	{
	case nullable_bool_t::True:
		o << L"True";
		break;
	case nullable_bool_t::False:
		o << L"False";
		break;
	case nullable_bool_t::Undefined:
		o << L"Undefined";
		break;
	default:
		o << L"unknown nullable bool value";
		break;
	}
	return o;
}


/// @brief Преобразует nullable_bool в строку
inline std::string nullable_bool_to_str(nullable_bool_t nb) {
    std::string str_value;
    if (nb == nullable_bool_t::True) {
        str_value = "True";
    }
    else if (nb == nullable_bool_t::False) {
        str_value = "False";
    }
    else {
        str_value = "Undefined";
    }
    return str_value;
}


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

