#pragma once


TEST(Fixed, LinearConstraints1) {

    typedef typename fixed_system_types<2>::var_type var_type;

    var_type x{ 0, 0 };
    var_type dx{ 1, 1 };

    fixed_linear_constraints<2, 1> linear_constraints;
    linear_constraints.a = { 1 , 1 };
    linear_constraints.b = 1;
    linear_constraints.trim(x, dx); // �������� � ������������

    // dx ������ ���������� � {1,1} �� {0.5, 0.5}

}

TEST(Fixed, LinearConstraints2) {

    typedef typename fixed_system_types<2>::var_type var_type;

    var_type x{ 0, 0 };
    var_type dx{ 0, 1 };

    fixed_linear_constraints<2, 1> linear_constraints;
    linear_constraints.a = { -1 , 1 };
    linear_constraints.b = 0;

    linear_constraints.trim(x, dx); // ���������� ����� �� �����������
    // dx ������ ��������������� � {0,1} �� {0.5, 0.5}
}

TEST(Fixed, LinearConstraints3) {

    typedef typename fixed_system_types<2>::var_type var_type;

    var_type x{ 0, 0 };
    var_type dx{ 0, 1 };

    fixed_linear_constraints<2, 1> linear_constraints;
    linear_constraints.a = { 1 , 1 };
    linear_constraints.b = 0;

    linear_constraints.trim(x, dx); // ���������� ����� �� �����������
    // dx ������ ��������������� � {0,1} �� {-0.5, 0.5}

    var_type ethalon{ -0.5, 0.5 };
}

TEST(Fixed, LinearConstraints4) {

    typedef typename fixed_system_types<2>::var_type var_type;

    var_type x{ 0, 0 };
    var_type dx{ 1, 0 };

    fixed_linear_constraints<2, 1> linear_constraints;
    linear_constraints.a = { 1 , 1 };
    linear_constraints.b = 0;

    linear_constraints.trim(x, dx); // ���������� ����� �� �����������

    var_type ethalon{ 0.5, -0.5 };

    // dx ������ ��������������� � {1, 0} �� {0.5, -0.5}
}

TEST(Fixed, LinearConstraints5) {

    typedef typename fixed_system_types<2>::var_type var_type;

    var_type x{ 0, 0 };
    var_type dx{ 1, 1 };

    fixed_linear_constraints<2, 1> linear_constraints;
    linear_constraints.a = { 1 , 1 };
    linear_constraints.b = 0;

    linear_constraints.trim(x, dx); // ���������� ����� �� �����������
    // dx ������ ��������������� � {1, 1} �� {0, 0}

    var_type ethalon{ 0.0, 0.0 };

}