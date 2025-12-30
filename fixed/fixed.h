#pragma once

#ifndef __FIXED_H__
#define __FIXED_H__

#define _USE_MATH_DEFINES // подключить константы
#include <cmath>
#include <numeric>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "fixed/helpers/math_helpers.h"
#include "fixed/helpers/string_helpers.h"

namespace fixed_solvers {
;

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
