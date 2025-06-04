#pragma once

/// �������� �������� ���������� �� ������������
inline constexpr double eps_constraints = 1e-8;

template <std::ptrdiff_t Dimension>
struct fixed_solver_constraints;

/// @brief �������� ������ ��� ������� �� ������������
/// @param callback ���������� � ����������� index, min, max. 
/// ���� ����������� ������
template <typename Callback>
inline void prepare_box_constraints(
    const std::vector<std::pair<size_t, double>>& minimum,
    const std::vector<std::pair<size_t, double>>& maximum,
    Callback callback
)
{
    size_t iimin = 0;
    size_t iimax = 0;

    bool have_min = iimin < minimum.size();
    bool have_max = iimax < maximum.size();
    while (have_min || have_max)
    {
        size_t min_index = std::numeric_limits<size_t>::max();
        size_t max_index = std::numeric_limits<size_t>::max();

        if (have_min) {
            min_index = minimum[iimin].first;
        }
        if (have_max) {
            max_index = maximum[iimax].first;
        }

        double min_value = -std::numeric_limits<double>::infinity();
        double max_value = std::numeric_limits<double>::infinity();

        size_t var_index;
        if (have_min && have_max) {
            if (min_index < max_index) {
                var_index = min_index;
                min_value = minimum[iimin++].second;
            }
            else if (max_index < min_index) {
                var_index = max_index;
                max_value = maximum[iimax++].second;
            }
            else {
                var_index = min_index;
                min_value = minimum[iimin++].second;
                max_value = maximum[iimax++].second;
            }
        }
        else if (have_min) {
            var_index = min_index;
            min_value = minimum[iimin++].second;
        }
        else { // have_max
            var_index = max_index;
            max_value = maximum[iimax++].second;
        }
        callback(var_index, min_value, max_value);
        //minqpsetbci(state, var_index, min_value, max_value);

        have_min = iimin < minimum.size();
        have_max = iimax < maximum.size();
    }

}


/*! \brief C�������� ��������� ����������� ��� �������� ���������� �����������
*
* ����� ����, ��� �� ����� �� ������������, �� �� ������ � Doxygen ������������ */
template <>
struct fixed_solver_constraints<-1>
{
    /// @brief ������ ����������� �� ������������� ���������� ���������
    std::vector<std::pair<size_t, double>> relative_boundary;
    /// @brief ������ ����������� �� ����������� �������� ���������
    std::vector<std::pair<size_t, double>> minimum;
    /// @brief ������ ����������� �� ������������ �������� ���������
    std::vector<std::pair<size_t, double>> maximum;
    /// @brief ���������� ���������� ����������� 
    size_t get_constraint_count() const
    {
        return minimum.size() + maximum.size();
    }

    std::pair<
        std::vector<std::pair<size_t, double>>,
        std::vector<std::pair<size_t, double>>
    > get_relative_constraints(
        const VectorXd& current_argument) const
    {
        // �����������
        auto n = current_argument.size();

        std::vector<std::pair<size_t, double>> mins = this->minimum;
        std::vector<std::pair<size_t, double>> maxs = this->maximum;

        for (auto& [index, min_value] : mins) {
            //x0 + dx > min   ==>   dx > min - x0;
            min_value -= current_argument(index);
        }

        for (auto& [index, max_value] : maxs) {
            //x0 + dx < max   ==>   dx < max - x0;
            max_value -= current_argument(index);
        }

        return std::make_pair(std::move(mins), std::move(maxs));
    }


    /// @brief ����������� �� �������� � �������� ��� ������������� ����������������
    static std::pair<MatrixXd, VectorXd> get_inequalities_constraints_vectors_dense(
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
    std::pair<MatrixXd, VectorXd> get_inequalities_constraints_dense(const size_t argument_size) const
    {
        // �����������
        const auto n = argument_size;
        MatrixXd A = MatrixXd::Zero(get_constraint_count(), n);
        VectorXd B = VectorXd::Zero(get_constraint_count());

        size_t offset = 0;
        // ������������ ��������
        {
            auto [a, b] = get_inequalities_constraints_vectors_dense(n, maximum);

            A.block(offset, 0, a.rows(), n) = a;
            B.segment(offset, b.size()) = b;
            offset += a.rows();
        }
        // ����������� ��������
        {
            auto [a, b] = get_inequalities_constraints_vectors_dense(n, minimum);

            A.block(offset, 0, a.rows(), n) = -a;
            B.segment(offset, b.size()) = -b;
            offset += a.rows();
        }
        return std::make_pair(A, B);
    }

    /// @brief ��������� �� ����������
    void trim_relative(VectorXd& increment) const
    {
        using std::max;
        double factor = 1;
        for (const auto& kvp : relative_boundary)
        {
            int sign = sgn(increment(kvp.first));
            if (sign * increment(kvp.first) > kvp.second) {
                double current_factor = sign * increment(kvp.first) / kvp.second;
                factor = max(factor, current_factor);
            }
        }
        if (factor > 1)
        {
            increment /= factor;
        }
    }
    /// @brief ��������� ������� ���������� ���������, ����������� �� ������������ min ��� max
    /// �� �� relative, �.�. �� ��� ����������� ����� ������� ������ 
    /// @param argument 
    /// @return 
    bool has_active_constraints(const VectorXd& argument) const {
        for (const auto& [index, min_value] : minimum) {
            if (std::abs(argument[index] - min_value) < eps_constraints) {
                return true;
            }
        }
        for (const auto& [index, max_value] : maximum) {
            if (std::abs(argument[index] - max_value) < eps_constraints) {
                return true;
            }
        }
        return false;
    }

public:
    std::pair<SparseMatrix<double, Eigen::ColMajor>, VectorXd> get_inequalities_constraints_sparse(
        const VectorXd& current_argument) const
    {
        // �����������
        auto n = current_argument.size();

        size_t row_index = 0;
        std::vector<Eigen::Triplet<double>> A;
        VectorXd b(minimum.size() + maximum.size());

        auto add_constraints = [&](const std::vector<std::pair<size_t, double>>& boundaries, double sign) {
            for (const auto& [var_index, value] : boundaries) {
                Eigen::Triplet<double> triplet(
                    static_cast<int>(row_index),
                    static_cast<int>(var_index),
                    sign * 1.0
                );

                A.emplace_back(std::move(triplet));
                b(row_index) = sign * value;
                row_index++;
            }
            };

        add_constraints(minimum, -1.0);
        add_constraints(maximum, +1.0);

        SparseMatrix<double, Eigen::ColMajor> A_matrix(minimum.size() + maximum.size(), n);
        A_matrix.setFromTriplets(A.begin(), A.end());
        b -= A_matrix * current_argument;

        return std::make_pair(std::move(A_matrix), std::move(b));
    }

    /// @brief ��������� �� ���������
    void trim_max(VectorXd& argument, VectorXd& increment) const
    {
        using std::max;
        double factor = 1;

        for (const auto& [index, max_value] : maximum)
        {
            if (argument[index] + increment[index] > max_value) {
                if (std::abs(argument[index] - max_value) < eps_constraints) {
                    // �������� ��� ��� �� �����������, allowed_decrement ����� �������,
                    // �������������� factor ���������� ����������� (��. ����� "else" ����)
                    // �� ��������� ��� ���������� ��� ������� factor, ����� �������� 
                    increment[index] = 0;
                }
                else {
                    double allowed_increment = max_value - argument[index];
                    factor = max(factor, abs(increment[index]) / allowed_increment);
                }
            }
        }
        if (factor > 1)
        {
            increment = increment / factor;
        }

    }

    /// @brief ��������� �� ��������
    void trim_min(VectorXd& argument, VectorXd& increment) const
    {
        using std::max;
        double factor = 1;

        for (const auto& [index, min_value] : minimum)
        {
            if (argument[index] + increment[index] < min_value) {
                if (std::abs(argument[index] - min_value) < eps_constraints) {
                    // �������� ��� ��� �� �����������, allowed_decrement ����� �������,
                    // �������������� factor ���������� ����������� (��. ����� "else" ����)
                    // �� ��������� ��� ���������� ��� ������� factor, ����� �������� 
                    increment[index] = 0;
                }
                else {
                    double allowed_decrement = argument[index] - min_value;
                    factor = max(factor, abs(increment[index]) / allowed_decrement);
                }

            }
        }
        if (factor > 1)
        {
            increment = increment / factor;
        }
    }

    /// @brief ��������� �� ������������
    void ensure_constraints(VectorXd& argument) const
    {
        VectorXd increment = VectorXd::Zero(argument.size());

        for (const auto& [index, min_value] : minimum) {
            if (argument[index] < min_value) {
                argument[index] = min_value;
            }
        }
        for (const auto& [index, max_value] : maximum) {
            if (argument[index] > max_value) {
                argument[index] = max_value;
            }
        }
    }
};



///@brief ������ ��������� ��������� ����������� ������� ������������� �����������
template <std::ptrdiff_t Dimension>
struct fixed_solver_constraints
{
    /// ��������� ������� ����������
    typedef typename fixed_system_types<Dimension>::var_type var_type;
    /// ����������� �� ����������
    var_type relative_boundary{ fixed_system_types<Dimension>::default_var() };
    /// ����������� �� ��������
    var_type minimum{ fixed_system_types<Dimension>::default_var() };
    /// ����������� �� ���������
    var_type maximum{ fixed_system_types<Dimension>::default_var() };


    std::pair<
        std::vector<std::pair<size_t, double>>,
        std::vector<std::pair<size_t, double>>
    > get_relative_constraints(
        const var_type& current_argument) const
    {
        throw std::runtime_error("not impl");
    }

    /*!
    * \brief ������� ���� �� �������� ��� ���������� ������
    *
    * \param [in] argument �������� ��������� �� ������� ��������
    * \param [in] increment �������� ���������� �� ������� ��������
    */
    void trim_min(double argument, double& increment) const
    {
        if (std::isnan(minimum))
            return;
        if (argument + increment < minimum) {
            increment = minimum - argument;
        }
    }

    /*!
    * \brief ������� ���� �� �������� ��� ���������� ������
    *
    * \param [in] argument �������� ��������� �� ������� ��������
    * \param [in] increment �������� ���������� �� ������� ��������
    */
    void trim_min(const array<double, Dimension>& argument,
        array<double, Dimension>& increment) const
    {
        using std::max;
        double factor = 1;

        for (size_t index = 0; index < increment.size(); ++index) {
            if (std::isnan(minimum[index])) {
                continue;
            }

            if (argument[index] + increment[index] < minimum[index]) {


                if (std::abs(argument[index] - minimum[index]) < eps_constraints) {
                    // �������� ��� ��� �� �����������, allowed_decrement ����� �������,
                    // �������������� factor ���������� ����������� (��. ����� "else" ����)
                    // �� ��������� ��� ���������� ��� ������� factor, ����� �������� 
                    increment[index] = 0;
                }
                else {
                    double allowed_decrement = argument[index] - minimum[index];
                    factor = max(factor, abs(increment[index]) / allowed_decrement);
                }

            }
        }
        if (factor > 1)
        {
            increment = increment / factor;
        }
    }

    /*!
    * \brief �������� �������� ��������� ������ ����������� ���/����
    *
    * \param [in] argument �������� ��������� �� ������� ��������
    */
    void ensure_constraints(double& argument) const
    {
        using std::min;
        using std::max;
        if (!std::isnan(maximum)) {
            argument = min(argument, maximum);
        }
        if (!std::isnan(minimum)) {
            argument = max(argument, minimum);
        }
    }
    void ensure_constraints(std::array<double, 2>& argument) const
    {
        using std::min;
        using std::max;
        for (size_t index = 0; index < Dimension; ++index) {
            if (!std::isnan(maximum[index])) {
                argument[index] = min(argument[index], maximum[index]);
            }
            if (!std::isnan(minimum[index])) {
                argument[index] = max(argument[index], minimum[index]);
            }

        }
    }

    bool has_active_constraints(const std::array<double, Dimension>& argument) const {
        for (size_t index = 0; index < Dimension; ++index) {
            if (!std::isnan(minimum[index])) {
                if (std::abs(argument[index] - minimum[index]) < eps_constraints) {
                    return true;
                }
            }
            if (!std::isnan(maximum[index])) {
                if (std::abs(argument[index] - maximum[index]) < eps_constraints) {
                    return true;
                }
            }
        }
        return false;
    }
    bool has_active_constraints(const double& argument) const {
        if (!std::isnan(minimum)) {
            if (std::abs(argument - minimum) < eps_constraints) {
                return true;
            }
        }
        if (!std::isnan(maximum)) {
            if (std::abs(argument - maximum) < eps_constraints) {
                return true;
            }
        }
        return false;
    }




    std::pair<SparseMatrix<double, Eigen::ColMajor>, VectorXd> get_inequalities_constraints_sparse(
        const std::array<double, Dimension>& current_argument) const
    {
        // �����������
        auto n = current_argument.size();


        std::vector<Eigen::Triplet<double>> A;
        std::vector<double> b;

        auto add_constraints = [&](const var_type& boundaries, double sign) {
            for (size_t var_index = 0; var_index < Dimension; ++var_index) {
                if (std::isnan(boundaries[var_index]))
                    continue;

                size_t row_index = A.size();

                Eigen::Triplet<double> triplet(
                    static_cast<int>(row_index),
                    static_cast<int>(var_index),
                    sign * 1.0
                );

                A.emplace_back(std::move(triplet));
                b.emplace_back(sign * boundaries[var_index]);
            }
            };

        add_constraints(minimum, -1.0);
        add_constraints(maximum, +1.0);

        SparseMatrix<double, Eigen::ColMajor> A_matrix(A.size(), n);
        A_matrix.setFromTriplets(A.begin(), A.end());

        VectorXd b_vec = VectorXd::Map(b.data(), b.size());
        b_vec -= A_matrix * VectorXd::Map(current_argument.data(), current_argument.size());

        return std::make_pair(std::move(A_matrix), std::move(b_vec));
    }

    /*!
    * \brief ������� ���� �� ��������� ��� ���������� ������
    *
    * \param [in] argument �������� ��������� �� ������� ��������
    * \param [in] increment �������� ���������� �� ������� ��������
    */
    void trim_max(double argument, double& increment) const
    {
        if (std::isnan(maximum))
            return;

        if (argument + increment > maximum) {
            increment = maximum - argument;
        }
    }

    /*!
    * \brief ������� ���� �� ��������� ��� ���������� ������
    *
    * \param [in] argument �������� ��������� �� ������� ��������
    * \param [in] increment �������� ���������� �� ������� ��������
    */
    void trim_max(const array<double, Dimension>& argument,
        array<double, Dimension>& increment) const
    {
        using std::max;
        double factor = 1;

        for (size_t index = 0; index < increment.size(); ++index) {
            if (std::isnan(maximum[index])) {
                continue;
            }

            if (argument[index] + increment[index] > maximum[index]) {
                if (std::abs(argument[index] - maximum[index]) < eps_constraints) {
                    // �������� ��� ��� �� �����������, allowed_increment ����� �������,
                    // �������������� factor ���������� ����������� (��. ����� "else" ����)
                    // �� ��������� ��� ���������� ��� ������� factor, ����� �������� 
                    increment[index] = 0;
                }
                else {
                    double allowed_increment = maximum[index] - argument[index];
                    factor = max(factor, increment[index] / allowed_increment);
                }
            }
        }
        if (factor > 1)
        {
            increment = increment / factor;
        }
    }

    /*!
    * \brief ������� �� ���������� ��� ���������� ������
    *
    * \param [in] increment �������� ���������� ��������� �� ������� ��������
    */
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

    /*!
    * \brief ������� �� ���������� ��� ���������� ������
    *
    * \param [in] increment �������� ���������� ��������� �� ������� ��������
    */
    void trim_relative(array<double, Dimension>& increment) const
    {
        using std::max;
        double factor = 1;
        for (size_t index = 0; index < increment.size(); ++index) {
            if (std::isnan(relative_boundary[index])) {
                continue;
            }

            double abs_inc = abs(increment[index]);
            if (abs_inc > relative_boundary[index]) {
                double current_factor = abs_inc / relative_boundary[index];
                factor = max(factor, current_factor);
            }
        }
        if (factor > 1)
        {
            increment = increment / factor;
        }
    }
};


///@brief �������� ����������� - ������������ �����������
template <std::ptrdiff_t Dimension, std::ptrdiff_t Count>
struct fixed_linear_constraints;


///@brief ������������� ���������� �������� ����������� 
template <std::ptrdiff_t Dimension>
struct fixed_linear_constraints<Dimension, 0> {

    /// ��������� ������� ����������
    typedef typename fixed_system_types<Dimension>::var_type var_type;

    /// ������� trim ������ �� ������, �� ������ ����� �������
    inline void trim(const var_type& argument, var_type& increment) const
    {
    }
};

///@brief ������������� ��� ������ ������� �������, ���� �����������
template <>
struct fixed_linear_constraints<2, 1> {
    /// @brief ��� ���������
    typedef typename fixed_system_types<2>::var_type var_type;
    /// @brief ��� ������� �������
    typedef typename fixed_system_types<2>::matrix_type matrix_type;

    /// ������������ a (����� ����� ax <= b)
    var_type a{ std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN() };

    /// ������������ b (������ ����� ax <= b)
    double b{ std::numeric_limits<double>::quiet_NaN() };

    /*!
    * \brief ��������, ��� ����������� x �� �������� �����������
    *
    * \param [in] x ������� �����������
    */
    bool check_constraint_satisfaction(const var_type& x) const {
        if (std::isfinite(b))
            return inner_prod(a, x) <= b;
        else
            return true;
    }

    /*!
    * \brief ��������, ��� ����������� x ��������� �� ������� �����������
    *
    * \param [in] x ������� �����������
    */
    bool check_constraint_border(const var_type& x) const {
        if (std::isfinite(b))
            return std::abs(inner_prod(a, x) - b) < eps_constraints;
        else
            return true;
    }

    /// @brief �� ������ ������������� ��������� ������ ������ ������������ ���������
    /// ax = b
    /// @param p1 ������ �����
    /// @param p2  ������ �����
    /// @return ���� ������, ������: (a, b)
    static std::pair<var_type, double> get_line_coeffs(const var_type& p1, const var_type& p2) {
        double x1 = p1[0];
        double y1 = p1[1];
        double x2 = p2[0];
        double y2 = p2[1];

        // ��� � ���� y = kx + b
        double k = (y2 - y1) / (x2 - x1);
        double b = y1 - k * x1;

        // -k*x + 1y = b

        return make_pair(var_type{ -k, 1.0 }, b);
    }

    /// @brief ����� �� �������� ������������ a'x <= b
    /// @param argument ������� ��������
    /// @param increment ����������
    void trim(const var_type& x, var_type& dx) const
    {
        if (!std::isfinite(b))
            return;

        const var_type& p1 = x;
        const var_type p2 = x + dx;

        if (check_constraint_satisfaction(p2))
            return;

        if (check_constraint_border(p1)) {
            // ��� ���-�� ������ ����� ������ �� �������
            double k = -a[0] / a[1];
            double alpha = atan(k);
            double p = sqrt(dx[0] * dx[0] + dx[1] * dx[1]);
            double beta = acos(dx[0] / p);

            double gamma = beta - alpha;
            double p_dash = p * cos(gamma);
            double px = p_dash * cos(alpha);
            double py = p_dash * sin(alpha);
            dx = { px, py };
        }
        else {
            auto [a2, b2] = get_line_coeffs(p1, p2);
            var_type x_star = solve_linear_system({ a, a2 }, { b, b2 });
            dx = x_star - x;
        }
    }
};

/// @brief �������� ��� �����������, �� ���� ��������� ������ �������, 
/// �.�. ����������� ����� trim
template <std::ptrdiff_t Dimension, std::ptrdiff_t Count>
struct fixed_linear_constraints
{
};
