#pragma once



TEST(Common, LinearEquationSolver)
{
    array<array<double, 2>, 2> A;
    A[0] = { 2, 2 };
    A[1] = { 3, 4 };

    array<double, 2> b = { 6, 11 };


    auto x = solve_linear_system(A, b);
}
