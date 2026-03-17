#pragma once

/// @brief Положительное значение меньше eps - должно быть заменено на +eps
TEST(EnsureAbsEpsilonValue, PositiveValueLessThanEpsilon) {
    // Arrange
    double input = 1e-10;
    // Act
    double result = fixed_solvers::ensure_abs_epsilon_value(input);
    // Assert
    double expected = 1e-6;
    EXPECT_DOUBLE_EQ(result, expected);
}

/// @brief Отрицательное значение меньше eps - должно быть заменено на -eps
TEST(EnsureAbsEpsilonValue, NegativeValueLessThanEpsilon) {
    // Arrange
    double input = -1e-10;
    // Act
    double result = fixed_solvers::ensure_abs_epsilon_value(input);
    // Assert
    double expected = -1e-6;
    EXPECT_DOUBLE_EQ(result, expected);
}

/// @brief Положительное значение больше eps - должно остаться без изменений
TEST(EnsureAbsEpsilonValue, PositiveValueGreaterThanEpsilon) {
    // Arrange
    double input = 1.0;
    // Act
    double result = fixed_solvers::ensure_abs_epsilon_value(input);
    // Assert
    double expected = 1.0;
    EXPECT_DOUBLE_EQ(result, expected);
}

/// @brief Отрицательное значение больше eps - должно остаться без изменений
TEST(EnsureAbsEpsilonValue, NegativeValueGreaterThanEpsilon) {
    // Arrange
    double input = -1.0;
    // Act
    double result = fixed_solvers::ensure_abs_epsilon_value(input);
    // Assert
    double expected = -1.0;
    EXPECT_DOUBLE_EQ(result, expected);
}
