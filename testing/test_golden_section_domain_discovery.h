#pragma once

namespace {
golden_section_domain_discovery_parameters make_domain_parameters(
    size_t iteration_count = 16,
    bool use_continuous_domain_based_proposal = false
) {
    golden_section_domain_discovery_parameters parameters;
    parameters.iteration_count = iteration_count;
    parameters.function_decrement_factor = std::numeric_limits<double>::quiet_NaN();
    parameters.function_target_value = std::numeric_limits<double>::quiet_NaN();
    parameters.use_continuous_domain_based_proposal = use_continuous_domain_based_proposal;
    return parameters;
}
}

/// @brief Возвращает минимум в стандартном режиме, когда все точки в области определения.
TEST(GoldenSectionDomainDiscovery, ReturnsMinimumWhenAllPointsInDomain)
{
    // Arrange: настраиваем параметры и квадратичную функцию с минимумом в 0.6.
    auto parameters = make_domain_parameters();
    auto function = [](double x) { return std::pow(x - 0.6, 2.0); };

    // Act: запускаем поиск минимума на [0, 1].
    auto [step, iterations] = golden_section_search_domain_discovery::search(
        parameters, function, 0.0, 1.0, function(0.0), function(1.0));

    // Assert: проверяем сходимость к ожидаемому минимуму.
    ASSERT_TRUE(std::isfinite(step));
    ASSERT_GT(iterations, 0u);
    ASSERT_NEAR(step, 0.6, 1e-2);
}

/// @brief Выбирает правую часть интервала, когда alpha вне области определения, а beta определена.
TEST(GoldenSectionDomainDiscovery, ChoosesRightPartWhenAlphaOutAndBetaIn)
{
    // Arrange: задаем функцию с разрывом ООФ и минимумом справа от разрыва.
    auto parameters = make_domain_parameters(14);
    auto function = [](double x) {
        const bool in_left = (x >= 0.0) && (x < 0.2);
        const bool in_right = (x > 0.6) && (x <= 1.0);
        if (!in_left && !in_right) {
            throw domain_violation{};
        }
        return std::pow(x - 0.9, 2.0);
    };

    // Act: запускаем поиск, при котором alpha попадает в разрыв ООФ.
    auto [step, iterations] = golden_section_search_domain_discovery::search(
        parameters, function, 0.0, 1.0, function(0.0), function(1.0));

    // Assert: проверяем смещение результата в правую часть интервала.
    ASSERT_TRUE(std::isfinite(step));
    ASSERT_GT(iterations, 0u);
    ASSERT_GT(step, 0.6);
}

/// @brief Возвращает fail-контракт, когда при alpha вне области определения минимум не выводится из унимодальности.
TEST(GoldenSectionDomainDiscovery, ReturnsFailContractWhenAlphaOutAndBetaInWithoutUnimodalityEvidence)
{
    // Arrange: задаем функцию, где значения по обе стороны разрыва ООФ не дают вывода из унимодальности.
    auto parameters = make_domain_parameters(14);
    auto function = [](double x) {
        const bool in_left = (x >= 0.0) && (x < 0.2);
        const bool in_right = (x > 0.6) && (x <= 1.0);
        if (!in_left && !in_right) {
            throw domain_violation{};
        }
        if (in_left) {
            return std::pow(x - 0.03, 2.0);
        }
        return 10.0 + std::pow(x - 0.9, 2.0);
    };

    // Act: запускаем поиск, ожидая исчерпание итераций без локализации.
    auto [step, iterations] = golden_section_search_domain_discovery::search(
        parameters, function, 0.0, 1.0, function(0.0), function(1.0));

    // Assert: проверяем fail-контракт по step=NaN и счетчику итераций.
    ASSERT_FALSE(std::isfinite(step));
    ASSERT_EQ(iterations, parameters.iteration_count + 1);
}

/// @brief Выбирает левую часть интервала, когда beta вне области определения и правая граница определена.
TEST(GoldenSectionDomainDiscovery, ChoosesLeftPartWhenBetaOutAndBDefined)
{
    // Arrange: задаем функцию с границей ООФ и явное f(b) для активации ветки 5-7.
    auto parameters = make_domain_parameters(18);
    const double domain_border = 0.5;
    auto function = [&](double x) {
        if (x >= domain_border) {
            throw domain_violation{};
        }
        return std::pow(x - 0.35, 2.0);
    };

    // Act: запускаем поиск с явным f(b), задающим defined_b == true.
    auto [step, iterations] = golden_section_search_domain_discovery::search(
        parameters, function, 0.0, 1.0, function(0.0), 0.01);

    // Assert: проверяем, что результат лежит в левой части и близок к минимуму 0.35.
    ASSERT_TRUE(std::isfinite(step));
    ASSERT_GT(iterations, 0u);
    ASSERT_LT(step, domain_border);
    ASSERT_NEAR(step, 0.35, 2e-2);
}

/// @brief Выбирает правую часть интервала, когда beta вне области определения и выполняется правило для сдвига a.
TEST(GoldenSectionDomainDiscovery, ChoosesRightPartWhenBetaOutAndBDefined)
{
    // Arrange: задаем функцию и f(b), при которых выполняется условие f(alpha) > f(b) (ветка 6).
    auto parameters = make_domain_parameters(10);
    auto function = [](double x) {
        if (x >= 0.5) {
            throw domain_violation{};
        }
        return std::pow(x - 0.45, 2.0);
    };

    // Act: запускаем поиск, ожидая сдвиг левой границы a := alpha.
    auto [step, iterations] = golden_section_search_domain_discovery::search(
        parameters, function, 0.0, 1.0, function(0.0), 1e-6);

    // Assert: проверяем, что найденный аргумент сместился вправо.
    ASSERT_TRUE(std::isfinite(step));
    ASSERT_GT(iterations, 0u);
    ASSERT_GT(step, 0.3);
}

/// @brief Бросает logic_error, когда beta вне области определения и минимум не принадлежит текущему интервалу.
TEST(GoldenSectionDomainDiscovery, ThrowsLogicErrorWhenBetaOutAndBDefinedAndMinimumNotInInterval)
{
    // Arrange: задаем строго монотонную функцию с f(b) > f(alpha) — минимум не внутри [a, b].
    auto parameters = make_domain_parameters(6);
    auto function = [](double x) {
        if (x >= 0.5) {
            throw domain_violation{};
        }
        return 0.1 * x;
    };

    // Act: запускаем поиск с монотонной функцией.
    // Assert: проверяем выброс std::logic_error.
    EXPECT_THROW(
        golden_section_search_domain_discovery::search(
            parameters, function, 0.0, 1.0, function(0.0), 0.1),
        std::logic_error
    );
}

/// @brief Использует эвристику ООФ, когда alpha и beta вне области определения.
TEST(GoldenSectionDomainDiscovery, UsesDomainHeuristicWhenAlphaAndBetaOut)
{
    // Arrange: задаем функцию с границей ООФ и включаем эвристику связной ООФ.
    auto parameters = make_domain_parameters(18, true);
    const double domain_border = 0.38;
    auto function = [&](double x) {
        if (x >= domain_border) {
            throw domain_violation{};
        }
        return std::pow(x - 0.15, 2.0);
    };

    // Act: запускаем поиск без явного f(b) для активации ветки эвристики.
    auto [step, iterations] = golden_section_search_domain_discovery::search(
        parameters, function, 0.0, 1.0, function(0.0), std::numeric_limits<double>::quiet_NaN());

    // Assert: проверяем, что результат локализован внутри ООФ.
    ASSERT_TRUE(std::isfinite(step));
    ASSERT_GT(iterations, 0u);
    ASSERT_LT(step, domain_border);
}

/// @brief Возвращает fail-контракт, когда beta и b вне области определения и эвристика отключена.
TEST(GoldenSectionDomainDiscovery, ReturnsFailContractWhenBetaOutAndBOutWithoutHeuristic)
{
    // Arrange: задаем функцию с границей ООФ и выключаем эвристику.
    auto parameters = make_domain_parameters(12, false);
    const double domain_border = 0.55;
    auto function = [&](double x) {
        if (x >= domain_border) {
            throw domain_violation{};
        }
        return std::pow(x - 0.2, 2.0);
    };

    // Act: запускаем поиск без явного f(b), ожидая исчерпание итераций.
    auto [step, iterations] = golden_section_search_domain_discovery::search(
        parameters, function, 0.0, 1.0, function(0.0), std::numeric_limits<double>::quiet_NaN());

    // Assert: проверяем fail-контракт по step=NaN и счетчику итераций.
    ASSERT_FALSE(std::isfinite(step));
    ASSERT_EQ(iterations, parameters.iteration_count + 1);
}

/// @brief Продолжает поиск, когда beta и b вне области определения и эвристика включена.
TEST(GoldenSectionDomainDiscovery, ContinuesSearchWhenBetaOutAndBOutWithHeuristic)
{
    // Arrange: задаем функцию с границей ООФ и включаем эвристику связной ООФ.
    auto parameters = make_domain_parameters(18, true);
    const double domain_border = 0.55;
    auto function = [&](double x) {
        if (x >= domain_border) {
            throw domain_violation{};
        }
        return std::pow(x - 0.2, 2.0);
    };

    // Act: запускаем поиск без явного f(b).
    auto [step, iterations] = golden_section_search_domain_discovery::search(
        parameters, function, 0.0, 1.0, function(0.0), std::numeric_limits<double>::quiet_NaN());

    // Assert: проверяем, что результат локализован внутри ООФ.
    ASSERT_TRUE(std::isfinite(step));
    ASSERT_GT(iterations, 0u);
    ASSERT_LT(step, domain_border);
}

/// @brief Возвращает fail-контракт, когда функция вернула NaN в процессе поиска.
TEST(GoldenSectionDomainDiscovery, ReturnsFailContractWhenFunctionReturnsNan)
{
    // Arrange: задаем функцию, возвращающую NaN на части интервала (не domain_violation).
    auto parameters = make_domain_parameters(10);
    auto function = [](double x) {
        if (x > 0.6) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return std::pow(x - 0.2, 2.0);
    };

    // Act: запускаем поиск, ожидая распространения NaN в результат.
    auto [step, iterations] = golden_section_search_domain_discovery::search(
        parameters, function, 0.0, 1.0, function(0.0), std::numeric_limits<double>::quiet_NaN());

    // Assert: проверяем fail-контракт по step=NaN и счетчику итераций.
    ASSERT_FALSE(std::isfinite(step));
    ASSERT_EQ(iterations, parameters.iteration_count + 1);
}
