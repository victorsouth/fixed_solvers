#pragma once

struct linear_search_parameters {
    /// @brief ������� ��������� <1, ����� �� ���� ������������
    double maximum_step{ 1 };
    double minimum_step{ 0.05 };
    double step_divider{ 1.5 };
};



/// @brief ����������� �������
template <std::ptrdiff_t Dimension>
struct fixed_solver_constraints;

/// @brief ����������� �������
template <>
struct fixed_solver_constraints<-1>
{
    std::vector<std::pair<size_t, double>> relative_boundary;
    std::vector<std::pair<size_t, double>> minimum;
    std::vector<std::pair<size_t, double>> maximum;
    size_t get_constraint_count() const
    {
        return minimum.size() + maximum.size();
    }

    static std::pair<MatrixXd, VectorXd> get_inequalities_constraints_vectors(
        size_t argument_dimension, 
        const std::vector<std::pair<size_t, double>>& boundaries)
    {
        MatrixXd A = MatrixXd::Zero(boundaries.size(), argument_dimension);
        VectorXd b = VectorXd::Zero(boundaries.size());

        int row_index = 0;
        for (const auto& kvp : boundaries) {
            A(row_index, kvp.first) = 1;
            b(row_index) = kvp.second;
            row_index++;
        }
        return std::make_pair(A, b);
    }

    /// @brief ��������� ������ minimum, maximum
    std::pair<MatrixXd, VectorXd> get_inequalities_constraints(const size_t argument_size) const
    {
        // �����������
        const auto n = argument_size;
        MatrixXd A = MatrixXd::Zero(get_constraint_count(), n);
        VectorXd B = VectorXd::Zero(get_constraint_count());

        size_t offset = 0;
        // ������������ ��������
        { 
            auto [a, b] = get_inequalities_constraints_vectors(n, maximum);

            A.block(offset, 0, a.rows(), n) = a;
            B.segment(offset, b.size()) = b;
            offset += a.rows();
        }
        // ����������� ��������
        {
            auto [a, b] = get_inequalities_constraints_vectors(n, minimum);

            A.block(offset, 0, a.rows(), n) = -a;
            B.segment(offset, b.size()) = -b;
            offset += a.rows();
        }
        return std::make_pair(A, B);
    }


    void trim_increment_relative(VectorXd& increment) const
    {
        double factor = 1;
        for (const auto& kvp : relative_boundary)
        {
            int sign = sgn(increment(kvp.first));
            if (sign * increment(kvp.first) > kvp.second) {
                double current_factor = sign * increment(kvp.first) / kvp.second;
                factor = std::max(factor, current_factor);
            }
        }
        if (factor > 1)
        {
            increment /= factor;
        }
    }

};


/// @brief ����������� �������
template <std::ptrdiff_t Dimension>
struct fixed_solver_constraints
{
    typedef typename fixed_system_types<Dimension>::var_type var_type;
    var_type relative_boundary{ fixed_system_types<Dimension>::default_var() };
    var_type minimum{ fixed_system_types<Dimension>::default_var() };
    var_type maximum{ fixed_system_types<Dimension>::default_var() };

    void trim_relative(double& increment) const
    {
        if (std::isnan(relative_boundary))
            return;
        double abs_inc = abs(increment);
        if (abs_inc > relative_boundary) {
            double factor = abs_inc / relative_boundary;
            increment /= factor;
        }
    }

    void trim_relative(array<double, Dimension>& increment) const
    {
        double factor = 1;
        for (size_t index = 0; index < increment.size(); ++index) {
            if (std::isnan(relative_boundary[index])) {
                continue;
            }

            double abs_inc = abs(increment[index]);
            if (abs_inc > relative_boundary[index]) {
                double current_factor = abs_inc / relative_boundary[index];
                factor = std::max(factor, current_factor);
            }
        }
        if (factor > 1)
        {
            increment = increment / factor;
        }
    }
};


/// @brief ��������� ��������� 
template <std::ptrdiff_t Dimension>
struct fixed_solver_parameters_t
{
    fixed_solver_constraints<Dimension> constraints;
    linear_search_parameters line_search;
    bool step_criteria_assuming_search_step{ false };
    double argument_increment_norm{ 1e-4 };
    size_t iteration_count{ 100 };
};

enum class numerical_result_code_t
{
    NoNumericalError, IllConditionedMatrix, LargeConditionNumber,
    NotConverged, NumericalNanValues, Converged
};

/// @brief ��������� ������� ���������� ������
template <std::ptrdiff_t Dimension>
struct fixed_solver_result_t {
    typedef typename fixed_system_types<Dimension>::var_type var_type;
    typedef typename fixed_system_types<Dimension>::right_party_type function_type;

    /// @brief ������� ����������
    double argument_increment_metric{ 0 };
    /// @brief ����������, ��� ���������� ����������� ���������� ������
    bool argument_increment_criteria{ false };
    /// @brief ���������� ������������ ���������
    numerical_result_code_t result_code{ numerical_result_code_t::NotConverged };
    /// @brief ���������� ������� �� ��������� ���������� ������
    function_type residuals;
    /// @brief ������� �������� �� ��������� ���������� ������
    var_type argument;
    /// @brief ���������� ��������
    size_t iteration_count{ 0 };
};


/// @brief ������ �������-������� ��� ������ ��������� ������������� �����������
/// � ��� ����� ���������
/// �������-������� - �������� ����������� ����
/// ����������� ����������� ���� �� ���� ���������
template <size_t Dimension>
class fixed_newton_raphson {
public:
    typedef typename fixed_system_types<Dimension>::var_type var_type;
    typedef typename fixed_system_types<Dimension>::right_party_type function_type;
    typedef typename fixed_system_types<Dimension>::equation_coeffs_type equation_coeffs_type;
public:
    static inline bool has_not_finite(const double value)
    {
        if (!std::isfinite(value)) {
            return true;
        }
        return false;
    }

    template <typename Container>
    static inline bool has_not_finite(const Container& values)
    {
        for (const double value : values) {
            if (has_not_finite(value))
                return true;
        }
        return false;
    }

public:
    template <typename Function>
    static inline double line_search(Function& function,
        const linear_search_parameters& parameters)
    {

        //double function_initial = function(alpha_prev);
        bool found_better_than_initial = false;

        double function_initial = function(0);

        double alpha_prev, alpha_curr;
        double function_prev, function_current;

        // ����� ������������ �������� alpha, ��� ������� ������ �� �������������
        alpha_curr = parameters.maximum_step;
        function_current = function(alpha_curr);

        while (alpha_curr >= parameters.minimum_step)
        {
            alpha_prev = alpha_curr;
            alpha_curr /= parameters.step_divider;
            function_prev = function_current;
            function_current = function(alpha_curr);

            if (function_current < function_initial)
                found_better_than_initial = true;

            if (found_better_than_initial) {
                // ����, ����� �.�. ������ �����
                if (function_current > function_prev) {
                    return alpha_prev;
                }
                else
                    continue;
            }
        }

        return std::numeric_limits<double>::quiet_NaN();
    }

    static double argument_increment_factor(
        const var_type& argument, const var_type& argument_increment);


    /// @brief ��������� ������� ������ ��������� ������� ������-�������
    /// @param residuals ������� �������
    /// @param initial_argument ��������� �����������
    /// @param solver_parameters ��������� �������
    /// @param result ���������� �������
    static void solve_dense(
        fixed_system_t<Dimension>& residuals,
        const var_type& initial_argument,
        const fixed_solver_parameters_t<Dimension>& solver_parameters,
        fixed_solver_result_t<Dimension>* result
    )
    {
        var_type& r = result->residuals;
        function_type& argument = result->argument;
        double& argument_increment_metric = result->argument_increment_metric;
        bool& argument_increment_criteria = result->argument_increment_criteria;
        

        argument = initial_argument;
        var_type argument_increment;
        var_type p;

        r = residuals.residuals(argument);
        if (has_not_finite(r)) {
            r = residuals.residuals(argument); // ��� �������
            result->result_code = numerical_result_code_t::NumericalNanValues;
            return;
        }

        result->result_code = numerical_result_code_t::NotConverged;

        size_t& iteration = result->iteration_count;
        for (iteration = 0; iteration < solver_parameters.iteration_count; ++iteration)
        {
            auto J = residuals.jacobian_dense(argument);

            p = -solve_linear_system(J, r); // ����� ������ ���� ���������� ������

            // todo: ���������� ������

            if (solver_parameters.step_criteria_assuming_search_step == false)
            {
                argument_increment_metric = argument_increment_factor(argument, p);
                argument_increment_criteria =
                    argument_increment_metric < solver_parameters.argument_increment_norm;
                if (argument_increment_criteria) {
                    result->result_code = numerical_result_code_t::Converged;
                    break;
                }
            }

            solver_parameters.constraints.trim_relative(p);

            auto directed_function = [&](double step) {
                return residuals(argument + step * p);
            };

            double search_step = line_search(directed_function, solver_parameters.line_search);
            if (std::isnan(search_step)) {
                result->result_code = numerical_result_code_t::NotConverged;
                break;
            }

            argument_increment = search_step * p;

            argument += argument_increment;

            r = residuals.residuals(argument);
            if (has_not_finite(r)) {
                r = residuals.residuals(argument); // ��� �������
                result->result_code = numerical_result_code_t::NumericalNanValues;
                break;
            }

            argument_increment_metric = solver_parameters.step_criteria_assuming_search_step
                ? argument_increment_factor(argument, argument_increment)
                : argument_increment_factor(argument, p);
            argument_increment_criteria =
                argument_increment_metric < solver_parameters.argument_increment_norm;
            bool custom_criteria = false; 
            if (custom_criteria || argument_increment_criteria) {
                result->result_code = numerical_result_code_t::Converged;
                break;
            }
        }
    }

};

template <>
inline double fixed_newton_raphson<1>::argument_increment_factor(
    const double& argument,
    const double& argument_increment)
{
    double arg = std::max(1.0, argument);
    double inc = argument_increment;

    double result = abs(inc / arg);
    return result;
}

template <size_t Dimension>
inline double fixed_newton_raphson<Dimension>::argument_increment_factor(
    const var_type& argument, const var_type& argument_increment)
{
    double squared_sum = 0;

    for (int component = 0; component < argument_increment.size(); ++component) {
        double arg = std::max(1.0, argument[component]);
        double inc = argument_increment[component];
        squared_sum += pow(inc / arg, 2);
    }
    return sqrt(squared_sum) / Dimension;
}