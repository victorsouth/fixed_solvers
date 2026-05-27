#pragma once

/// @brief Возвращает alpha=0 без итераций, когда обе граничные ц.ф. ниже function_target_value.
TEST(GoldenSectionTargetValue, ReturnsZeroStep_WhenBothBoundaryValuesBelowTarget)
{
    // Arrange: задаем порог и функцию с малыми значениями на концах отрезка.
    golden_section_parameters parameters;
    parameters.iteration_count = 10;
    parameters.function_decrement_factor = std::numeric_limits<double>::quiet_NaN();
    parameters.function_target_value = 1e-8;
    const double f_a = 1e-10;
    const double f_b = 1e-11;
    auto function = [](double x) {
        return std::pow(x - 0.5, 2.0);
    };

    // Act: запускаем поиск с известными f(a) и f(b).
    auto [step, iterations] = golden_section_search::search(
        parameters, function, 0.0, 1.0, f_a, f_b);

    // Assert: шаг Ньютона не делается, итерации ЗС не выполняются.
    ASSERT_EQ(step, 0.0);
    ASSERT_EQ(iterations, 0);
}

/// @brief Возвращает меньшую границу без итераций, когда только одна граничная ц.ф. ниже порога.
TEST(GoldenSectionTargetValue, ReturnsSmallerBoundary_WhenOnlyOneBoundaryValueBelowTarget)
{
    // Arrange: f(a) ниже порога, f(b) выше порога, минимум на a.
    golden_section_parameters parameters;
    parameters.iteration_count = 10;
    parameters.function_decrement_factor = std::numeric_limits<double>::quiet_NaN();
    parameters.function_target_value = 1e-8;
    const double f_a = 1e-10;
    const double f_b = 1.0;
    auto function = [](double x) {
        return std::pow(x - 0.5, 2.0);
    };

    // Act: запускаем поиск с известными f(a) и f(b).
    auto [step, iterations] = golden_section_search::search(
        parameters, function, 0.0, 1.0, f_a, f_b);

    // Assert: возвращается левая граница без итераций ЗС.
    ASSERT_EQ(step, 0.0);
    ASSERT_EQ(iterations, 0);
}

/// @brief Возвращает правую границу без итераций, когда только f(b) ниже порога.
TEST(GoldenSectionTargetValue, ReturnsRightBoundary_WhenOnlyRightBoundaryValueBelowTarget)
{
    // Arrange: f(b) ниже порога, f(a) выше порога.
    golden_section_parameters parameters;
    parameters.iteration_count = 10;
    parameters.function_decrement_factor = std::numeric_limits<double>::quiet_NaN();
    parameters.function_target_value = 1e-8;
    const double f_a = 1.0;
    const double f_b = 1e-10;
    auto function = [](double x) {
        return std::pow(x - 0.5, 2.0);
    };

    // Act: запускаем поиск с известными f(a) и f(b).
    auto [step, iterations] = golden_section_search::search(
        parameters, function, 0.0, 1.0, f_a, f_b);

    // Assert: возвращается правая граница без итераций ЗС.
    ASSERT_EQ(step, 1.0);
    ASSERT_EQ(iterations, 0);
}

/// @brief Запускает итерации ЗС, когда обе граничные ц.ф. выше function_target_value.
TEST(GoldenSectionTargetValue, RunsIterations_WhenBothBoundaryValuesAboveTarget)
{
    // Arrange: обе граничные ц.ф. выше порога, минимум внутри отрезка.
    golden_section_parameters parameters;
    parameters.iteration_count = 3;
    parameters.function_decrement_factor = std::numeric_limits<double>::quiet_NaN();
    parameters.function_target_value = 1e-8;
    auto function = [](double x) {
        return std::pow(x - 0.3, 2.0) + 0.1;
    };

    // Act: запускаем поиск.
    auto [step, iterations] = golden_section_search::search(
        parameters, function, 0.0, 1.0, function(0.0));

    // Assert: поиск выполняет итерации и находит конечный шаг.
    ASSERT_TRUE(std::isfinite(step));
    ASSERT_GT(iterations, 0);
}

/// @brief Возвращает alpha=0 без итераций в domain_discovery при обеих граничных ц.ф. ниже порога.
TEST(GoldenSectionTargetValue, DomainDiscoveryReturnsZeroStep_WhenBothBoundaryValuesBelowTarget)
{
    // Arrange: порог и малые значения на концах отрезка.
    golden_section_domain_discovery_parameters parameters;
    parameters.iteration_count = 10;
    parameters.function_decrement_factor = std::numeric_limits<double>::quiet_NaN();
    parameters.function_target_value = 1e-8;
    parameters.mode = domain_discovery_mode_t::require_connected_domain;
    const double f_a = 1e-10;
    const double f_b = 1e-11;
    auto function = [](double x) {
        return std::pow(x - 0.5, 2.0);
    };

    // Act: запускаем поиск с известными f(a) и f(b).
    auto [step, iterations] = golden_section_search_domain_discovery::search(
        parameters, function, 0.0, 1.0, f_a, f_b);

    // Assert: шаг Ньютона не делается, итерации ЗС не выполняются.
    ASSERT_EQ(step, 0.0);
    ASSERT_EQ(iterations, 0);
}
