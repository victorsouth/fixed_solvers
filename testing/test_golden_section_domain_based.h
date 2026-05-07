#pragma once

/// @brief Проверяет штатный сценарий золотого сечения, когда все рабочие точки лежат в ООФ.
TEST(GoldenSectionDomainBased, HandlesCase1AllPointsInDomain)
{
    // Arrange: настраиваем параметры и унимодальную функцию на полном интервале.
    golden_section_parameters parameters;
    parameters.iteration_count = 16;
    parameters.function_decrement_factor = std::numeric_limits<double>::quiet_NaN();
    parameters.function_target_value = std::numeric_limits<double>::quiet_NaN();

    auto function = [](double x) {
        return std::pow(x - 0.6, 2.0);
    };

    double a = 0.0;
    double b = 1.0;

    // Act: запускаем domain-based поиск.
    auto [step, iterations] = golden_section_search::search(
        parameters, function, a, b, function(a), function(b));

    // Assert: шаг конечен, итерации выполнены, найдено положение минимума.
    ASSERT_TRUE(std::isfinite(step));
    ASSERT_GT(iterations, 0u);
    ASSERT_NEAR(step, 0.6, 1e-2);
}

/// @brief Проверяет сценарий, когда правая граница интервала вне ООФ, но внутренние точки доступны.
TEST(GoldenSectionDomainBased, HandlesCase2WhenRightBorderOutOfDomain)
{
    // Arrange: задаем правую границу ООФ и функцию с domain_violation за границей.
    golden_section_parameters parameters;
    parameters.iteration_count = 16;
    parameters.function_decrement_factor = std::numeric_limits<double>::quiet_NaN();
    parameters.function_target_value = std::numeric_limits<double>::quiet_NaN();

    const double domain_border = 0.92;
    auto function = [&](double x) {
        if (x >= domain_border) {
            throw domain_violation{};
        }
        return std::pow(x - 0.85, 2.0);
    };

    double a = 0.0;
    double b = 1.0;

    // Act: запускаем поиск при неизвестном f(b).
    auto [step, iterations] = golden_section_search::search(
        parameters, function, a, b, function(a), std::numeric_limits<double>::quiet_NaN());

    // Assert: найден корректный конечный шаг внутри ООФ.
    ASSERT_TRUE(std::isfinite(step));
    ASSERT_GT(iterations, 0u);
    ASSERT_LT(step, domain_border);
}

/// @brief Проверяет сценарий, когда beta выходит за ООФ и поиск сужается влево.
TEST(GoldenSectionDomainBased, HandlesCase3WhenBetaOutOfDomain)
{
    // Arrange: задаем связную ООФ с границей внутри [a, b].
    golden_section_parameters parameters;
    parameters.iteration_count = 18;
    parameters.function_decrement_factor = std::numeric_limits<double>::quiet_NaN();
    parameters.function_target_value = std::numeric_limits<double>::quiet_NaN();

    const double domain_border = 0.5;
    auto function = [&](double x) {
        if (x >= domain_border) {
            throw domain_violation{};
        }
        return std::pow(x - 0.35, 2.0);
    };

    double a = 0.0;
    double b = 1.0;

    // Act: выполняем поиск при недоступной правой части интервала.
    auto [step, iterations] = golden_section_search::search(
        parameters, function, a, b, function(a), std::numeric_limits<double>::quiet_NaN());

    // Assert: шаг остается в ООФ и близок к ожидаемому минимуму.
    ASSERT_TRUE(std::isfinite(step));
    ASSERT_GT(iterations, 0u);
    ASSERT_LT(step, domain_border);
    ASSERT_NEAR(step, 0.35, 2e-2);
}

/// @brief Проверяет случай, когда при разрывной ООФ выбирается правый интервал [beta, b].
TEST(GoldenSectionDomainBased, HandlesCase4AsDisconnectedDomainFailure)
{
    // Arrange: формируем несвязную ООФ и делаем f(beta) <= f(a) для перехода вправо.
    golden_section_parameters parameters;
    parameters.iteration_count = 12;
    parameters.function_decrement_factor = std::numeric_limits<double>::quiet_NaN();
    parameters.function_target_value = std::numeric_limits<double>::quiet_NaN();

    auto function = [](double x) {
        bool in_left = (x >= 0.0 && x < 0.2);
        bool in_right = (x > 0.6 && x <= 1.0);
        if (!in_left && !in_right) {
            throw domain_violation{};
        }
        return std::pow(x - 0.9, 2.0);
    };

    double a = 0.0;
    double b = 1.0;

    // Act: запускаем поиск на интервале с разрывом ООФ.
    auto [step, iterations] = golden_section_search::search(
        parameters, function, a, b, function(a), function(b));

    // Assert: поиск выбирает правую связную часть и дает конечный шаг.
    ASSERT_TRUE(std::isfinite(step));
    ASSERT_GT(iterations, 0u);
    ASSERT_GT(step, 0.6);
}

/// @brief Проверяет случай, когда при разрывной ООФ выбирается левый интервал [a, alpha].
TEST(GoldenSectionDomainBased, HandlesCase4WhenLeftIntervalChosen)
{
    // Arrange: формируем несвязную ООФ и делаем f(beta) > f(a) для перехода влево.
    golden_section_parameters parameters;
    parameters.iteration_count = 12;
    parameters.function_decrement_factor = std::numeric_limits<double>::quiet_NaN();
    parameters.function_target_value = std::numeric_limits<double>::quiet_NaN();

    auto function = [](double x) {
        bool in_left = (x >= 0.0 && x < 0.2);
        bool in_right = (x > 0.6 && x <= 1.0);
        if (!in_left && !in_right) {
            throw domain_violation{};
        }
        if (in_left) {
            return std::pow(x - 0.05, 2.0);
        }
        return 10.0 + std::pow(x - 0.9, 2.0);
    };

    double a = 0.0;
    double b = 1.0;

    // Act: запускаем поиск на интервале с разрывом ООФ.
    auto [step, iterations] = golden_section_search::search(
        parameters, function, a, b, function(a), function(b));

    // Assert: поиск выбирает левую связную часть и дает конечный шаг.
    ASSERT_TRUE(std::isfinite(step));
    ASSERT_GT(iterations, 0u);
    ASSERT_LT(step, 0.2);
}

/// @brief Проверяет сценарий, когда обе внутренние точки alpha и beta вне ООФ.
TEST(GoldenSectionDomainBased, HandlesCase5WhenAlphaAndBetaOutOfDomain)
{
    // Arrange: задаем границу ООФ так, чтобы рабочие точки выпадали правее границы.
    golden_section_parameters parameters;
    parameters.iteration_count = 18;
    parameters.function_decrement_factor = std::numeric_limits<double>::quiet_NaN();
    parameters.function_target_value = std::numeric_limits<double>::quiet_NaN();

    const double domain_border = 0.38;
    auto function = [&](double x) {
        if (x >= domain_border) {
            throw domain_violation{};
        }
        return std::pow(x - 0.15, 2.0);
    };

    double a = 0.0;
    double b = 1.0;

    // Act: запускаем поиск и даем алгоритму сузиться к левой доступной части.
    auto [step, iterations] = golden_section_search::search(
        parameters, function, a, b, function(a), std::numeric_limits<double>::quiet_NaN());

    // Assert: получаем конечный шаг внутри ООФ.
    ASSERT_TRUE(std::isfinite(step));
    ASSERT_GT(iterations, 0u);
    ASSERT_LT(step, domain_border);
}

/// @brief Проверяет, что NaN в значении функции трактуется как численная ошибка поиска.
TEST(GoldenSectionDomainBased, TreatsFunctionNanAsSearchFailure)
{
    // Arrange: задаем функцию, возвращающую NaN на части интервала.
    golden_section_parameters parameters;
    parameters.iteration_count = 10;
    parameters.function_decrement_factor = std::numeric_limits<double>::quiet_NaN();
    parameters.function_target_value = std::numeric_limits<double>::quiet_NaN();

    auto function = [](double x) {
        if (x > 0.6) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return std::pow(x - 0.2, 2.0);
    };

    double a = 0.0;
    double b = 1.0;

    // Act: запускаем поиск при неизвестном f(b).
    auto [step, iterations] = golden_section_search::search(
        parameters, function, a, b, function(a), std::numeric_limits<double>::quiet_NaN());

    // Assert: получаем fail-контракт поиска.
    ASSERT_FALSE(std::isfinite(step));
    ASSERT_EQ(iterations, parameters.iteration_count + 1);
}

/// @brief Проверяет, что для неунимодальной функции выбрасывается исключение логики.
TEST(GoldenSectionDomainBased, ThrowsOnNonUnimodalFunction)
{
    // Arrange: задаем функцию с локальным максимумом в опорных точках ЗС.
    golden_section_parameters parameters;
    parameters.iteration_count = 8;
    parameters.function_decrement_factor = std::numeric_limits<double>::quiet_NaN();
    parameters.function_target_value = std::numeric_limits<double>::quiet_NaN();

    auto non_unimodal = [](double x) {
        // Локальный максимум около alpha на первом шаге ЗС.
        return -std::pow(x - 0.38, 2.0);
    };

    double a = 0.0;
    double b = 1.0;

    // Act: выполняем поиск на неунимодальной функции.
    // Assert: ожидаем исключение о нарушении предположения унимодальности.
    EXPECT_THROW(
        golden_section_search::search(
            parameters, non_unimodal, a, b, non_unimodal(a), non_unimodal(b)),
        std::logic_error
    );
}
