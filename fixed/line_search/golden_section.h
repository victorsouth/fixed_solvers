#pragma once

/// @brief ��������� ������ �������� �������
struct golden_section_parameters {
    /// @brief ������������ ��� �������
    double maximum_step{ 1 };
    /// @brief ���������� ��������
    size_t iteration_count{ 10 };
    /// @brief ��������� ���������� �������� �.�. � ���������
    double function_decrement_factor{ 2 };
    /// @brief ����������� (���������� �����) �������� ������� �������� �.�.
    /// ���� NaN, �� �� ������������
    double function_target_value{ std::numeric_limits<double>::quiet_NaN() };
    /// @brief ���������� ������������� ����� ������� ����������� ��������
    double get_final_section_length() const {
        return pow(0.618, iteration_count);
    }
    /// @brief ��������� ������� ���������� �� �������� �.�.
    bool decrement_factor_criteria(double f_current_min, double f_0) const {
        if (std::isfinite(function_decrement_factor)) {
            double decrement = f_0 / f_current_min;
            // ���������� ������� ����������
            return  decrement > function_decrement_factor;
        }
        else
            return false;
    }
    /// @brief ��������� ����������� �������� �.�.
    bool target_value_criteria(double f_current_min) const {
        return std::isfinite(function_target_value) && (f_current_min < function_target_value);
    }
};


/// @brief �������� ����� ������� �������� �������
class golden_section_search {
public:
    typedef golden_section_parameters parameters_type;
public:
    static inline double get_alpha(double a, double b) {
        return a + 2 / (3 + sqrt(5)) * (b - a);
    }
    static inline double get_beta(double a, double b) {
        return a + 2 / (1 + sqrt(5)) * (b - a);
    }
public:
    template <typename Function>
    static inline std::pair<double, size_t> search(
        const golden_section_parameters& parameters,
        Function& function,
        double a, double b,
        double f_a, double f_b
    )
    {
        auto check_convergence = [&](double f_min, double f_0) {
            return parameters.decrement_factor_criteria(f_min, f_0)
                || parameters.target_value_criteria(f_min);
        };

        // ��������� �������� �.�., �� �������� ������� ������ ��������
        double f_0 = f_a;

        // �������������� ������� ������
        auto [x_min, f_min] = f_b < f_a
            ? make_pair(b, f_b)
            : make_pair(a, f_a);
        // ���� �������� ����� ����������, �������
        if (check_convergence(f_min, f_0)) {
            // ��������, ��� ������ �� ������������ ��������
            return { x_min, 0 };
        }

        // ������� ��������, ���� ���������� ��� �������� f(alpha), f(beta)
        double alpha = get_alpha(a, b);
        double f_alpha = function(alpha);

        double beta = get_beta(a, b);
        double f_beta = function(beta);

        // ������� ��������� ������ x_min, f_min �� ��������� ��������� a, alpha, beta, b
        auto update_minimum = [&]() {
            // ����� �������� ������� �� ��������
            if (f_alpha < f_beta) {
                // 1. ������ �� ��������������, ������� ����� � ��������� [a, beta]
                // 2. �� ����� ��������� ���������� �������� f(a), f(alpha), f(beta)
                // ����� ����� ������� �� ���
                // 3. ��� ���� ��� ���������, ��� f(alpha) < f(beta)
                // �������� ����� min(f(alpha), f(a))
                std::tie(x_min, f_min) = f_alpha < f_a
                    ? make_pair(alpha, f_alpha)
                    : make_pair(a, f_a);
            }
            else {
                // 1. ������� � ��������� [alpha, b]
                // 2. �� ����� ��������� ���������� �������� f(alpha), f(beta), f(b)
                // 3. ��� ���� ��� ���������, ��� f(alpha) > f(beta)
                // �������� ����� min(f(b), f(beta))
                std::tie(x_min, f_min) = f_b < f_beta
                    ? make_pair(b, f_b)
                    : make_pair(beta, f_beta);
            }
        };
        update_minimum();
        if (check_convergence(f_min, f_0)) {
            return { x_min, 1 };
        }

        // ���� �� n-1 ��������, �.�. ������� �������� ���������� �� ��������� � ��� �������
        for (size_t index = 1; index < parameters.iteration_count; ++index) {
            if (f_alpha < f_beta) {
                // ��������� � [a, beta]
                // ������� [a, alpha] - "�������" ������� �������� �������
                b = beta; // b ��������, a ��������
                f_b = f_beta;
                beta = alpha;
                f_beta = f_alpha;
                alpha = get_alpha(a, b);
                f_alpha = function(alpha);
            }
            else {
                // ��������� � [alpha, b]. 
                // ������� [beta, b] - "�������" ������� �������� �������
                a = alpha; // a ��������, b ��������
                f_a = f_alpha;
                alpha = beta;
                f_alpha = f_beta;
                beta = get_beta(a, b);
                f_beta = function(beta);
            }
            update_minimum();
            if (check_convergence(f_min, f_0)) {
                return { x_min, index + 1 };
            }
        }

        if (f_min < f_0) {
            // �� ������� �� ��������� ��������, �� ������ ���������� ����
            // ����� ���������� �� ������� ���������� �� ��������� ��������, ������ n + 1
            return { x_min, parameters.iteration_count + 1 };
        }
        else {
            return {
                std::numeric_limits<double>::quiet_NaN(),
                parameters.iteration_count + 1
            };
        }
    }
};

