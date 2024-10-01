#pragma once


TEST(Fixed, LinearConstraints1) {

    typedef typename fixed_system_types<2>::var_type var_type;

    var_type x{ 0, 0 };
    var_type dx{ 1, 1 };

    fixed_linear_constraints<2, 1> linear_constraints;
    linear_constraints.a = { 1 , 1 };
    linear_constraints.b = 1;
    linear_constraints.trim(x, dx); // подходим к ограничениям

    // dx должно обрезаться с {1,1} до {0.5, 0.5}

}

TEST(Fixed, LinearConstraints2) {

    typedef typename fixed_system_types<2>::var_type var_type;

    var_type x{ 0, 0 };
    var_type dx{ 0, 1 };

    fixed_linear_constraints<2, 1> linear_constraints;
    linear_constraints.a = { -1 , 1 };
    linear_constraints.b = 0;

    linear_constraints.trim(x, dx); // изначально сидим на ограниченях
    // dx должно спроецироваться с {0,1} до {0.5, 0.5}
}

TEST(Fixed, LinearConstraints3) {

    typedef typename fixed_system_types<2>::var_type var_type;

    var_type x{ 0, 0 };
    var_type dx{ 0, 1 };

    fixed_linear_constraints<2, 1> linear_constraints;
    linear_constraints.a = { 1 , 1 };
    linear_constraints.b = 0;

    linear_constraints.trim(x, dx); // изначально сидим на ограниченях
    // dx должно спроецироваться с {0,1} до {-0.5, 0.5}

    var_type ethalon{ -0.5, 0.5 };
}

TEST(Fixed, LinearConstraints4) {

    typedef typename fixed_system_types<2>::var_type var_type;

    var_type x{ 0, 0 };
    var_type dx{ 1, 0 };

    fixed_linear_constraints<2, 1> linear_constraints;
    linear_constraints.a = { 1 , 1 };
    linear_constraints.b = 0;

    linear_constraints.trim(x, dx); // изначально сидим на ограниченях

    var_type ethalon{ 0.5, -0.5 };

    // dx должно спроецироваться с {1, 0} до {0.5, -0.5}
}

TEST(Fixed, LinearConstraints5) {

    typedef typename fixed_system_types<2>::var_type var_type;

    var_type x{ 0, 0 };
    var_type dx{ 1, 1 };

    fixed_linear_constraints<2, 1> linear_constraints;
    linear_constraints.a = { 1 , 1 };
    linear_constraints.b = 0;

    linear_constraints.trim(x, dx); // изначально сидим на ограниченях
    // dx должно спроецироваться с {1, 1} до {0, 0}

    var_type ethalon{ 0.0, 0.0 };

}