#pragma once

#ifndef __FIXED_H__
#define __FIXED_H__

#ifndef __COMMON_CUDA__

/// @brief Возвращает знак числа
/// @return для отрицательных -1, для положительных +1
template <typename T> inline int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

/// @brief Возведение вещественной переменной в квадрат
inline double sqr(double x) {
    return x * x;
}
#endif

#include <iomanip>
#include "fixed/qp/qp_wrapper.h"
#include "fixed/array_ext.h"
#include "fixed/fixed_system.h"
#include "fixed/fixed_linear_solver.h"
#include "fixed/fixed_nonlinear_solver.h"
#include "fixed/fixed_optimizer.h"

#endif
