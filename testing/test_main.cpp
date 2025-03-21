﻿
#define GTEST_BREAK_ON_FAILURE 1
#define GTEST_CATCH_EXCEPTIONS 0
#define GTEST_HAS_SEH 0
#define _VARIADIC_MAX 10 /* for gtest */
#include "gtest/gtest.h"

#include <fixed/fixed.h>

#include "test_linear_constraints.h"
#include "test_linear_system.h"
#include "test_newton_raphson.h"

//#include <iostream>


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
#ifdef _WIN32
    std::wcout.imbue(std::locale("rus_rus.866"));
#endif
int res= RUN_ALL_TESTS();
    return res;
}
