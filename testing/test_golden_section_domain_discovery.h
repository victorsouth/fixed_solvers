#pragma once

/// @brief Возвращает fail-контракт, когда при alpha вне ООФ и beta в ООФ интервал не удается локализовать.
TEST(GoldenSectionDomainDiscovery, ReturnsFailContract_WhenAlphaOutOfDomain_AndBetaInDomain_AndIntervalCannotBeLocalized)
{
    // Arrange: задаем функцию, где значения по обе стороны разрыва ООФ не дают вывода из унимодальности.
    golden_section_domain_discovery_parameters parameters;
    parameters.iteration_count = 14;
    parameters.function_decrement_factor = std::numeric_limits<double>::quiet_NaN();
    parameters.function_target_value = std::numeric_limits<double>::quiet_NaN();
    parameters.mode = domain_discovery_mode_t::require_connected_domain;
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
        parameters, function, 0.0, 1.0, function(0.0));

    // Assert: проверяем fail-контракт по step=NaN и счетчику итераций.
    ASSERT_FALSE(std::isfinite(step));
    ASSERT_EQ(iterations, parameters.iteration_count + 1);
}

/// @brief Бросает logic_error, когда beta и b вне области определения и требуется связность ООФ.
TEST(GoldenSectionDomainDiscovery, ThrowsLogicError_WhenDisconnectedDomainDetected)
{
    // Arrange: задаем функцию с разрывом ООФ и настаиваем на связной ООФ.
    golden_section_domain_discovery_parameters parameters;
    parameters.iteration_count = 12;
    parameters.function_decrement_factor = std::numeric_limits<double>::quiet_NaN();
    parameters.function_target_value = std::numeric_limits<double>::quiet_NaN();
    parameters.mode = domain_discovery_mode_t::require_connected_domain;
    const double domain_border = 0.55;
    auto function = [&](double x) {
        if (x >= domain_border) {
            throw domain_violation{};
        }
        return std::pow(x - 0.2, 2.0);
    };

    // Act: запускаем поиск в режиме только связной ООФ.
    // Assert: при обнаружении несвязной ООФ ожидаем logic_error.
    EXPECT_THROW(
        golden_section_search_domain_discovery::search(
            parameters, function, 0.0, 1.0, function(0.0)),
        std::logic_error
    );
}

/// @brief Бросает runtime_error, когда включен режим forbid_exit и функция нарушает ООФ.
TEST(GoldenSectionDomainDiscovery, ThrowsRuntimeError_WhenModeIsForbidExit_AndDomainViolationOccurs)
{
    // Arrange: включаем запрет выхода за ООФ и делаем b вне ООФ.
    golden_section_domain_discovery_parameters parameters;
    parameters.iteration_count = 8;
    parameters.function_decrement_factor = std::numeric_limits<double>::quiet_NaN();
    parameters.function_target_value = std::numeric_limits<double>::quiet_NaN();
    parameters.mode = domain_discovery_mode_t::forbid_exit;

    const double domain_border = 0.5;
    auto function = [&](double x) {
        if (x >= domain_border) {
            throw domain_violation{};
        }
        return std::pow(x - 0.1, 2.0);
    };

    // Act: запускаем поиск в режиме forbid_exit.
    // Assert: если b оценивается вне ООФ, ожидаем runtime_error.
    EXPECT_THROW(
        golden_section_search_domain_discovery::search(
            parameters, function, 0.0, 1.0, function(0.0)),
        std::runtime_error
    );
}

/// @brief Возвращает fail-контракт, когда функция вернула NaN в процессе поиска.
TEST(GoldenSectionDomainDiscovery, ReturnsFailContract_WhenFunctionReturnsNan)
{
    // Arrange: задаем функцию, возвращающую NaN на части интервала (не domain_violation).
    golden_section_domain_discovery_parameters parameters;
    parameters.iteration_count = 10;
    parameters.function_decrement_factor = std::numeric_limits<double>::quiet_NaN();
    parameters.function_target_value = std::numeric_limits<double>::quiet_NaN();
    parameters.mode = domain_discovery_mode_t::require_connected_domain;
    auto function = [](double x) {
        if (x > 0.6) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return std::pow(x - 0.2, 2.0);
    };

    // Act: запускаем поиск, ожидая распространения NaN в результат.
    auto [step, iterations] = golden_section_search_domain_discovery::search(
        parameters, function, 0.0, 1.0, function(0.0));

    // Assert: проверяем fail-контракт по step=NaN и счетчику итераций.
    ASSERT_FALSE(std::isfinite(step));
    ASSERT_EQ(iterations, parameters.iteration_count + 1);
}

/// @brief Возвращает конечный шаг для пункта 1: все точки в ООФ и f(alpha) < f(beta).
TEST(GoldenSectionDomainDiscovery, ReturnsFiniteStep_WhenAllInDomain_AndAlphaLessThanBeta_Rule01)
{
    // Arrange: все точки в ООФ, минимум в левой части.
    golden_section_domain_discovery_parameters parameters;
    parameters.iteration_count = 1;
    parameters.function_decrement_factor = std::numeric_limits<double>::quiet_NaN();
    parameters.function_target_value = std::numeric_limits<double>::quiet_NaN();
    parameters.mode = domain_discovery_mode_t::require_connected_domain;
    auto function = [](double x) { return std::pow(x - 0.2, 2.0); };

    // Act: запускаем один шаг поиска.
    auto [step, iterations] = golden_section_search_domain_discovery::search(
        parameters, function, 0.0, 1.0, function(0.0));

    // Assert: шаг найден и лежит в левой части.
    ASSERT_TRUE(std::isfinite(step));
    ASSERT_EQ(iterations, parameters.iteration_count + 1);
    ASSERT_LT(step, 0.7);
}

/// @brief Возвращает конечный шаг для пункта 2: все точки в ООФ и f(alpha) > f(beta).
TEST(GoldenSectionDomainDiscovery, ReturnsFiniteStep_WhenAllInDomain_AndAlphaGreaterThanBeta_Rule02)
{
    // Arrange: все точки в ООФ, минимум в правой части.
    golden_section_domain_discovery_parameters parameters;
    parameters.iteration_count = 1;
    parameters.function_decrement_factor = std::numeric_limits<double>::quiet_NaN();
    parameters.function_target_value = std::numeric_limits<double>::quiet_NaN();
    parameters.mode = domain_discovery_mode_t::require_connected_domain;
    auto function = [](double x) { return std::pow(x - 0.9, 2.0); };

    // Act: запускаем один шаг поиска.
    auto [step, iterations] = golden_section_search_domain_discovery::search(
        parameters, function, 0.0, 1.0, function(0.0));

    // Assert: шаг найден и смещен вправо.
    ASSERT_TRUE(std::isfinite(step));
    ASSERT_EQ(iterations, parameters.iteration_count + 1);
    ASSERT_GT(step, 0.3);
}

/// @brief Возвращает конечный шаг для пункта 3: b вне ООФ и f(alpha) < f(beta).
TEST(GoldenSectionDomainDiscovery, ReturnsFiniteStep_WhenBOutOfDomain_AndAlphaLessThanBeta_Rule03)
{
    // Arrange: b вне ООФ, alpha и beta в ООФ, минимум слева.
    golden_section_domain_discovery_parameters parameters;
    parameters.iteration_count = 1;
    parameters.function_decrement_factor = std::numeric_limits<double>::quiet_NaN();
    parameters.function_target_value = std::numeric_limits<double>::quiet_NaN();
    parameters.mode = domain_discovery_mode_t::require_connected_domain;
    auto function = [](double x) {
        if (x >= 0.99) {
            throw domain_violation{};
        }
        return std::pow(x - 0.2, 2.0);
    };

    // Act: запускаем один шаг без известного f(b).
    auto [step, iterations] = golden_section_search_domain_discovery::search(
        parameters, function, 0.0, 1.0, function(0.0));

    // Assert: шаг найден.
    ASSERT_TRUE(std::isfinite(step));
    ASSERT_EQ(iterations, parameters.iteration_count + 1);
}

/// @brief Возвращает конечный шаг для пункта 4: b вне ООФ и f(alpha) > f(beta).
TEST(GoldenSectionDomainDiscovery, ReturnsFiniteStep_WhenBOutOfDomain_AndAlphaGreaterThanBeta_Rule04)
{
    // Arrange: b вне ООФ, alpha и beta в ООФ, минимум справа.
    golden_section_domain_discovery_parameters parameters;
    parameters.iteration_count = 1;
    parameters.function_decrement_factor = std::numeric_limits<double>::quiet_NaN();
    parameters.function_target_value = std::numeric_limits<double>::quiet_NaN();
    parameters.mode = domain_discovery_mode_t::require_connected_domain;
    auto function = [](double x) {
        if (x >= 0.99) {
            throw domain_violation{};
        }
        return std::pow(x - 0.85, 2.0);
    };

    // Act: запускаем один шаг без известного f(b).
    auto [step, iterations] = golden_section_search_domain_discovery::search(
        parameters, function, 0.0, 1.0, function(0.0));

    // Assert: шаг найден.
    ASSERT_TRUE(std::isfinite(step));
    ASSERT_EQ(iterations, parameters.iteration_count + 1);
}

/// @brief Возвращает конечный шаг для пункта 5: beta вне ООФ и выбирается левая часть.
TEST(GoldenSectionDomainDiscovery, ReturnsFiniteStep_WhenBetaOutOfDomain_AndSelectsLeftPart_Rule05)
{
    // Arrange: разрыв ООФ между левой и правой компонентами; beta в разрыве, b в правой; f(alpha)<f(a).
    golden_section_domain_discovery_parameters parameters;
    parameters.iteration_count = 1;
    parameters.function_decrement_factor = std::numeric_limits<double>::quiet_NaN();
    parameters.function_target_value = std::numeric_limits<double>::quiet_NaN();
    parameters.mode = domain_discovery_mode_t::require_connected_domain;
    const double search_a = 0.0;
    const double search_b = 0.56;
    auto function = [](double x) {
        const bool in_left = (x >= 0.0) && (x < 0.22);
        const bool in_right = (x > 0.55) && (x <= 1.0);
        if (!in_left && !in_right) {
            throw domain_violation{};
        }
        if (in_left) {
            return std::pow(x - 0.18, 2.0);
        }
        return std::pow(x - 0.9, 2.0);
    };

    // Act: один шаг поиска; f(b) считается внутри search.
    auto [step, iterations] = golden_section_search_domain_discovery::search(
        parameters, function, search_a, search_b, function(search_a));

    // Assert: шаг найден.
    ASSERT_TRUE(std::isfinite(step));
    ASSERT_EQ(iterations, parameters.iteration_count + 1);
}

/// @brief Возвращает конечный шаг для пункта 6: beta вне ООФ и выбирается правая часть.
TEST(GoldenSectionDomainDiscovery, ReturnsFiniteStep_WhenBetaOutOfDomain_AndSelectsRightPart_Rule06)
{
    // Arrange: та же геометрия отрезка, что в Rule05; f(alpha)>f(b) — ветвь случая 6.
    golden_section_domain_discovery_parameters parameters;
    parameters.iteration_count = 1;
    parameters.function_decrement_factor = std::numeric_limits<double>::quiet_NaN();
    parameters.function_target_value = std::numeric_limits<double>::quiet_NaN();
    parameters.mode = domain_discovery_mode_t::allow_disconnected_domain;
    const double search_a = 0.0;
    const double search_b = 0.56;
    auto function = [](double x) {
        const bool in_left = (x >= 0.0) && (x < 0.22);
        const bool in_right = (x > 0.55) && (x <= 1.0);
        if (!in_left && !in_right) {
            throw domain_violation{};
        }
        if (in_left) {
            return 0.1 * x + 0.2;
        }
        return 1e-6;
    };

    // Act: один шаг поиска.
    auto [step, iterations] = golden_section_search_domain_discovery::search(
        parameters, function, search_a, search_b, function(search_a));

    // Assert: шаг найден.
    ASSERT_TRUE(std::isfinite(step));
    ASSERT_EQ(iterations, parameters.iteration_count + 1);
}

/// @brief Бросает logic_error для пункта 7: beta вне ООФ и минимум не в интервале.
TEST(GoldenSectionDomainDiscovery, ThrowsLogicError_WhenBetaOutOfDomain_AndMinimumNotInInterval_Rule07)
{
    // Arrange: та же геометрия; на левой компоненте монотонный рост, на правой низкий минимум — ветвь случая 7.
    golden_section_domain_discovery_parameters parameters;
    parameters.iteration_count = 1;
    parameters.function_decrement_factor = std::numeric_limits<double>::quiet_NaN();
    parameters.function_target_value = std::numeric_limits<double>::quiet_NaN();
    parameters.mode = domain_discovery_mode_t::require_connected_domain;
    const double search_a = 0.0;
    const double search_b = 0.56;
    auto function = [](double x) {
        const bool in_left = (x >= 0.0) && (x < 0.22);
        const bool in_right = (x > 0.55) && (x <= 1.0);
        if (!in_left && !in_right) {
            throw domain_violation{};
        }
        if (in_left) {
            return 0.1 * x;
        }
        return 1.0 + std::pow(x - 0.9, 2.0);
    };

    // Act: выполняем один шаг поиска.
    // Assert: ожидаем logic_error.
    EXPECT_THROW(
        golden_section_search_domain_discovery::search(
            parameters, function, search_a, search_b, function(search_a)),
        std::logic_error
    );
}

/// @brief Возвращает конечный шаг для пункта 8: alpha вне ООФ и выполняется переход в правую часть.
TEST(GoldenSectionDomainDiscovery, ReturnsFiniteStep_WhenAlphaOutOfDomain_AndSelectsRightPart_Rule08)
{
    // Arrange: alpha вне ООФ, beta и b в ООФ.
    golden_section_domain_discovery_parameters parameters;
    parameters.iteration_count = 1;
    parameters.function_decrement_factor = std::numeric_limits<double>::quiet_NaN();
    parameters.function_target_value = std::numeric_limits<double>::quiet_NaN();
    parameters.mode = domain_discovery_mode_t::allow_disconnected_domain;
    auto function = [](double x) {
        const bool in_left = (x >= 0.0) && (x < 0.2);
        const bool in_right = (x > 0.6) && (x <= 1.0);
        if (!in_left && !in_right) {
            throw domain_violation{};
        }
        return std::pow(x - 0.9, 2.0);
    };

    // Act: выполняем один шаг поиска.
    auto [step, iterations] = golden_section_search_domain_discovery::search(
        parameters, function, 0.0, 1.0, function(0.0));

    // Assert: шаг найден.
    ASSERT_TRUE(std::isfinite(step));
    ASSERT_EQ(iterations, parameters.iteration_count + 1);
}

/// @brief Возвращает конечный шаг для пункта 9: alpha вне ООФ и выбирается [a, beta].
TEST(GoldenSectionDomainDiscovery, ReturnsFiniteStep_WhenAlphaOutOfDomain_AndSelectsABeta_Rule09)
{
    // Arrange: alpha вне ООФ, beta в ООФ.
    golden_section_domain_discovery_parameters parameters;
    parameters.iteration_count = 1;
    parameters.function_decrement_factor = std::numeric_limits<double>::quiet_NaN();
    parameters.function_target_value = std::numeric_limits<double>::quiet_NaN();
    parameters.mode = domain_discovery_mode_t::allow_disconnected_domain;
    auto function = [](double x) {
        const bool in_left = (x >= 0.0) && (x < 0.2);
        const bool in_right = (x > 0.6) && (x <= 1.0);
        if (!in_left && !in_right) {
            throw domain_violation{};
        }
        if (in_left) {
            return 2.0 + std::pow(x - 0.05, 2.0);
        }
        return 0.1 + std::pow(x - 0.9, 2.0);
    };

    // Act: выполняем один шаг поиска.
    auto [step, iterations] = golden_section_search_domain_discovery::search(
        parameters, function, 0.0, 1.0, function(0.0));

    // Assert: шаг найден.
    ASSERT_TRUE(std::isfinite(step));
    ASSERT_EQ(iterations, parameters.iteration_count + 1);
}

/// @brief Бросает logic_error для пункта 10: alpha вне ООФ и минимум не в интервале.
TEST(GoldenSectionDomainDiscovery, ThrowsLogicError_WhenAlphaOutOfDomain_AndMinimumNotInInterval_Rule10)
{
    // Arrange: alpha вне ООФ, beta в ООФ и условие ошибки.
    golden_section_domain_discovery_parameters parameters;
    parameters.iteration_count = 1;
    parameters.function_decrement_factor = std::numeric_limits<double>::quiet_NaN();
    parameters.function_target_value = std::numeric_limits<double>::quiet_NaN();
    parameters.mode = domain_discovery_mode_t::allow_disconnected_domain;
    auto function = [](double x) {
        const bool in_left = (x >= 0.0) && (x < 0.2);
        const bool in_right = (x > 0.6) && (x <= 1.0);
        if (!in_left && !in_right) {
            throw domain_violation{};
        }
        if (in_left) {
            return 0.01 * x;
        }
        return 0.5 + std::pow(x - 0.6, 2.0);
    };

    // Act: выполняем один шаг поиска.
    // Assert: ожидаем logic_error.
    EXPECT_THROW(
        golden_section_search_domain_discovery::search(
            parameters, function, 0.0, 1.0, function(0.0)),
        std::logic_error
    );
}

/// @brief Возвращает конечный шаг для пункта 11: beta и b вне ООФ, выбирается [a, alpha].
TEST(GoldenSectionDomainDiscovery, ReturnsFiniteStep_WhenBetaAndBOutOfDomain_AndSelectsAAlpha_Rule11)
{
    // Arrange: beta и b вне ООФ, alpha в ООФ, эвристика включена.
    golden_section_domain_discovery_parameters parameters;
    parameters.iteration_count = 1;
    parameters.function_decrement_factor = std::numeric_limits<double>::quiet_NaN();
    parameters.function_target_value = std::numeric_limits<double>::quiet_NaN();
    parameters.mode = domain_discovery_mode_t::allow_disconnected_domain;
    auto function = [](double x) {
        if (x >= 0.55) {
            throw domain_violation{};
        }
        return std::pow(x - 0.2, 2.0);
    };

    // Act: выполняем один шаг поиска.
    auto [step, iterations] = golden_section_search_domain_discovery::search(
        parameters, function, 0.0, 1.0, function(0.0));

    // Assert: шаг найден.
    ASSERT_TRUE(std::isfinite(step));
    ASSERT_EQ(iterations, parameters.iteration_count + 1);
}

/// @brief Бросает logic_error для пункта 12 при выключенной эвристике.
TEST(GoldenSectionDomainDiscovery, ThrowsLogicError_WhenHeuristicDisabled_Rule12)
{
    // Arrange: beta и b вне ООФ, эвристика выключена.
    golden_section_domain_discovery_parameters parameters;
    parameters.iteration_count = 1;
    parameters.function_decrement_factor = std::numeric_limits<double>::quiet_NaN();
    parameters.function_target_value = std::numeric_limits<double>::quiet_NaN();
    parameters.mode = domain_discovery_mode_t::require_connected_domain;
    auto function = [](double x) {
        if (x >= 0.55) {
            throw domain_violation{};
        }
        return std::pow(x - 0.2, 2.0);
    };

    // Act: выполняем один шаг поиска.
    // Assert: в режиме связной ООФ при повторном выходе из ООФ ожидаем logic_error.
    EXPECT_THROW(
        golden_section_search_domain_discovery::search(
            parameters, function, 0.0, 1.0, function(0.0)),
        std::logic_error
    );
}

/// @brief Возвращает конечный шаг для пункта 13: alpha и b вне ООФ, выполняется [a, beta].
TEST(GoldenSectionDomainDiscovery, ReturnsFiniteStep_WhenAlphaAndBOutOfDomain_AndSelectsABeta_Rule13)
{
    // Arrange: alpha вне ООФ, beta в ООФ, b вне ООФ.
    golden_section_domain_discovery_parameters parameters;
    parameters.iteration_count = 1;
    parameters.function_decrement_factor = std::numeric_limits<double>::quiet_NaN();
    parameters.function_target_value = std::numeric_limits<double>::quiet_NaN();
    parameters.mode = domain_discovery_mode_t::allow_disconnected_domain;
    auto function = [](double x) {
        const bool in_left = (x >= 0.0) && (x < 0.2);
        const bool in_right = (x > 0.6) && (x < 0.9);
        if (!in_left && !in_right) {
            throw domain_violation{};
        }
        if (in_left) {
            return 2.0 + std::pow(x - 0.05, 2.0);
        }
        return 0.1 + std::pow(x - 0.8, 2.0);
    };

    // Act: выполняем один шаг поиска.
    auto [step, iterations] = golden_section_search_domain_discovery::search(
        parameters, function, 0.0, 1.0, function(0.0));

    // Assert: шаг найден.
    ASSERT_TRUE(std::isfinite(step));
    ASSERT_EQ(iterations, parameters.iteration_count + 1);
}

/// @brief Возвращает конечный шаг для пункта 14: alpha и b вне ООФ, применяется эвристика.
TEST(GoldenSectionDomainDiscovery, ReturnsFiniteStep_WhenAlphaAndBOutOfDomain_AndHeuristicEnabled_Rule14)
{
    // Arrange: alpha и b вне ООФ, beta в ООФ.
    golden_section_domain_discovery_parameters parameters;
    parameters.iteration_count = 1;
    parameters.function_decrement_factor = std::numeric_limits<double>::quiet_NaN();
    parameters.function_target_value = std::numeric_limits<double>::quiet_NaN();
    parameters.mode = domain_discovery_mode_t::allow_disconnected_domain;
    auto function = [](double x) {
        const bool in_left = (x >= 0.0) && (x < 0.2);
        const bool in_mid = (x > 0.55) && (x < 0.8);
        if (!in_left && !in_mid) {
            throw domain_violation{};
        }
        if (in_left) {
            return 2.0 + std::pow(x - 0.05, 2.0);
        }
        return 0.2 + std::pow(x - 0.62, 2.0);
    };

    // Act: выполняем один шаг поиска.
    auto [step, iterations] = golden_section_search_domain_discovery::search(
        parameters, function, 0.0, 1.0, function(0.0));

    // Assert: шаг найден за счет эвристики.
    ASSERT_TRUE(std::isfinite(step));
    ASSERT_EQ(iterations, parameters.iteration_count + 1);
}

/// @brief Возвращает конечный шаг для пункта 15: alpha и beta вне ООФ, используется эвристика.
TEST(GoldenSectionDomainDiscovery, ReturnsFiniteStep_WhenAlphaAndBetaOutOfDomain_AndHeuristicEnabled_Rule15)
{
    // Arrange: alpha и beta вне ООФ, a и b в ООФ.
    golden_section_domain_discovery_parameters parameters;
    parameters.iteration_count = 2;
    parameters.function_decrement_factor = std::numeric_limits<double>::quiet_NaN();
    parameters.function_target_value = std::numeric_limits<double>::quiet_NaN();
    parameters.mode = domain_discovery_mode_t::allow_disconnected_domain;
    auto function = [](double x) {
        const bool in_left = (x >= 0.0) && (x < 0.2);
        const bool in_right = (x > 0.95) && (x <= 1.0);
        if (!in_left && !in_right) {
            throw domain_violation{};
        }
        if (in_left) {
            return std::pow(x - 0.15, 2.0);
        }
        return 10.0 + std::pow(x - 1.0, 2.0);
    };

    // Act: выполняем один шаг поиска.
    auto [step, iterations] = golden_section_search_domain_discovery::search(
        parameters, function, 0.0, 1.0, function(0.0));

    // Assert: шаг найден за счет эвристики.
    ASSERT_TRUE(std::isfinite(step));
    ASSERT_EQ(iterations, parameters.iteration_count + 1);
}

/// @brief Возвращает конечный шаг для пункта 16: alpha, beta и b вне ООФ, используется эвристика.
TEST(GoldenSectionDomainDiscovery, ReturnsFiniteStep_WhenAlphaBetaAndBOutOfDomain_AndHeuristicEnabled_Rule16)
{
    // Arrange: alpha, beta и b вне ООФ, в ООФ только левая часть.
    golden_section_domain_discovery_parameters parameters;
    parameters.iteration_count = 3;
    parameters.function_decrement_factor = std::numeric_limits<double>::quiet_NaN();
    parameters.function_target_value = std::numeric_limits<double>::quiet_NaN();
    parameters.mode = domain_discovery_mode_t::allow_disconnected_domain;
    auto function = [](double x) {
        if (x >= 0.3) {
            throw domain_violation{};
        }
        return std::pow(x - 0.15, 2.0);
    };

    // Act: выполняем один шаг поиска.
    auto [step, iterations] = golden_section_search_domain_discovery::search(
        parameters, function, 0.0, 1.0, function(0.0));

    // Assert: шаг найден за счет эвристики.
    ASSERT_TRUE(std::isfinite(step));
    ASSERT_EQ(iterations, parameters.iteration_count + 1);
}
