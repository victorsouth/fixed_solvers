#pragma once
#include <fixed/fixed.h>
#include <iostream>

/// \brief метод решения уравнения
enum class fixed_bisectional_solution_type{
    /// метод дихотомии
    Bisection = 0 ,
    /// метод секущих
    Secant = 1,
    /// комбинированный настраиваемый в зависимомти от
    /// параметров решателя метод
    Combined = 2
};

/// \brief структура для хранения параметров решателя
struct fixed_bisectional_parameters_t{
    /// признак сохранения истории схождения
    bool argument_history{false};
    /// Собирает историю значений аргумента
    bool residual_history{ false };
    /// желаемая точность по аргументу
    double argument_precision{std::numeric_limits<double>::epsilon()};
    /// желаемая точность по невязке
    double residual_precision{std::numeric_limits<float>::epsilon()};
    /// тип рещения
    fixed_bisectional_solution_type solution_type = fixed_bisectional_solution_type::Bisection;
    /// нижняя граница итераций для метода секущих
    int secant_treshhold_iterations{3};
    /// нижняя граница по интервалу аргумента
    /// delta X min при которой считаем методом секущих
    /// если NaN. то её нет
    double secant_treshhold_min{std::numeric_limits<double>::quiet_NaN()};
    /// нижняя граница по интервалу аргумента
    /// delta X max при которой считаем методом секущих
    /// если NaN. то её нет
    double secant_treshhold_max{std::numeric_limits<double>::quiet_NaN()};
    /// минимальное значение аргумента -
    /// входной обязательный параметр
    double argument_limit_min{std::numeric_limits<double>::quiet_NaN()};
    /// максимальное значение аргумента -
    /// входной обязательный параметр
    double argument_limit_max{std::numeric_limits<double>::quiet_NaN()};
    /// признак печать отладочной информации
    /// в стандартный поток вывода
    bool verbose{false};
    /// использовать Illinois algorithm, сужает диапазон бисекции, ускоряя сходимость
    /// https://en.wikipedia.org/wiki/Regula_falsi
    bool use_Illinois{true};
    /// проверка схождения невязок на границе
    /// отключается в теплообменнике из-за непредсказуемости границ
    /// фактичеки поиск в интервале () вместо []
    bool check_boundary_before{true};
    /// @brief Проверка нахождения интервала поиска в диапазоне, 
    /// где разрешен метод секущих (secant_treshhold_max, secant_treshhold_min)
    bool secant_thresholds_satisfied(double x_lower, double x_upper) const 
    {
        double delta_current = std::abs(x_upper - x_lower);

        bool max_threshold_satisfied =
            std::isfinite(secant_treshhold_max) ? delta_current < secant_treshhold_max : true;
        bool min_threshold_satisfied =
            std::isfinite(secant_treshhold_min) ? delta_current > secant_treshhold_min : true;

        return max_threshold_satisfied && min_threshold_satisfied;
    }

};

/// @brief Результат расчета численного метода
template <std::ptrdiff_t Dimension>
struct fixed_bisection_result_t {
    /// \brief Псевдоним переменной
    /// определяем псевдоним искомой переменной для дальнейшего
    /// использования
    typedef typename fixed_system_types<Dimension>::var_type var_type;
    /// \brief Псевдоним невязки
    /// определяем Псевдоним невязки переменной для дальнейшего
    /// использования
    typedef typename fixed_system_types<Dimension>::right_party_type function_type;

    /// @brief Завершение итерационной процедуры
    numerical_result_code_t result_code{ numerical_result_code_t::NotConverged };
    /// @brief Балл сходимости
    convergence_score_t score{convergence_score_t::Poor};
    /// @brief Остаточная невязка по окончании численного метода
    function_type residuals{ fixed_system_types<Dimension>::default_var() };
    /// @brief Норма остаточной невязки по окончании численного метода
    double residuals_norm{0};
    /// @brief Искомый аргумент по окончании численного метода
    var_type argument{ fixed_system_types<Dimension>::default_var()};
    /// @brief Количество выполненных итераций
    size_t iteration_count{ 0 };
    /// @brief макс Количество выполненных итераций
    size_t max_allowed_iterations{ 0 };
    /// \brief достигнутая точность
    var_type reached_precision=std::numeric_limits<double>::epsilon();
};

/// @brief Аналитика сходимости
template <size_t Dimension>
struct fixed_bisection_result_analysis_t {
public:
    /// \brief Псевдоним переменной
    /// определяем псевдоним искомой переменной для дальнейшего
    /// использования
    typedef typename fixed_system_types<Dimension>::var_type var_type;
public:
    /// @brief Результат расчета
    std::vector<var_type> residual_history;
    /// @brief Результат расчета
    std::vector<var_type> argument_history;
};

/// \brief Солвер бисекции, а также прицепом к нему комбинированный солвер бисекции/секущих и отдельно солвер секущих
/// Логику выбора способа расчета шага см. в next_argument_value
template <size_t Dimension>
class fixed_bisectional {
public:
    /// \brief Псевдоним искомой переменной
    /// определяем Псевдоним искомой переменной для дальнейшего
    /// использования
    typedef typename fixed_system_types<Dimension>::var_type var_type;
    /// \brief Псевдоним невязки
    /// определяем Псевдоним невязки переменной для дальнейшего
    /// использования
    typedef typename fixed_system_types<Dimension>::right_party_type function_type;
private:
    /// @brief Возвращает количество итераций, нужных чтобы при делении пополам 
    /// изначальный отрезок initial_delta уменьшился до величины требуемой точности argument_precision
    /// Отрезок длинной initial_delta можно делить пополам 2^get_max_allowed_iterations раз, 
    /// пока он не станет равен argument_precision
    static size_t get_max_allowed_iterations(double initial_delta, double argument_precision) {
        return static_cast<size_t>(floor(log(1. / argument_precision) / log(2.0)));
    }

    /// @brief Проверяет, что диапазон поиска [x1, x2] меньше машинной точности x1 или x2
    /// Другими словами диапазон поиска меньше машинного эпсилон x1 или x2
    /// Эпсилон в IEE-754 это норма не для абсолютной погрешности, а для относительной!
    /// @param x1 Нижняя граница диапазона поиска
    /// @param x2 Верхняя граница диапазона поиска
    /// @return Истина, если диапазон настолько узок, что попал в машинный эпсилон - дальше улучшить не выйдет
    static bool is_within_machine_epsilon(double x1, double x2) {
        return std::abs(x1 - x2) <= std::max(std::abs(x2), std::abs(x1)) * std::numeric_limits<double>::epsilon();
    }

    /// \brief выдаёт следующее приближение аргумента в зависимости от типа расчета solver_parameters.solution_type
    /// \param solver_parameters - параметры солвера
    /// \param x1 - минимальное значение аргумента на предыдущем шаге
    /// \param x2 - минимальное значение аргумента на предыдущем шаге
    /// \param y1 - значение невязки при x1
    /// \param y2 - значение невязки при x2
    /// \param iterations - количество сделанных итераций
    /// \return возврвщяет следующее приближение аргумента, код ошибки расчета
    static std::pair<var_type, numerical_result_code_t> next_argument_value(const fixed_bisectional_parameters_t& solver_parameters,
        const var_type& x1,const var_type& x2,const var_type& y1,const var_type& y2,size_t iterations,bool& use_secant)
    {
        if(solver_parameters.verbose)
            std::cout<<"next_argument_value d:"<< std::abs(x2 - x1) <<std::endl;
        use_secant = false;

        switch (solver_parameters.solution_type)
        { 
        case fixed_bisectional_solution_type::Bisection:
            return next_argument_bisection(x1, x2,y1,y2);
        case fixed_bisectional_solution_type::Secant: {
            use_secant = true;
            return next_argument_secant(x1, x2, y1, y2);
        }
        case fixed_bisectional_solution_type::Combined: {
            use_secant = solver_parameters.secant_thresholds_satisfied(x1, x2)
                && static_cast<int>(iterations) > solver_parameters.secant_treshhold_iterations;
            if (use_secant) 
            {
                auto [x3, code] = next_argument_secant(x1, x2, y1, y2);
                if (code != numerical_result_code_t::Converged)  {
                    use_secant = false;
                    return next_argument_bisection(x1, x2,y1,y2);
                }
                else {
                    return std::make_pair(x3, code);
                }
            }
            else {
                return next_argument_bisection(x1, x2,y1,y2);
            }
        }
        default:
            throw std::runtime_error("unsupported solution_type");
        }
    }
    /// @brief Расчет следующего приближения по методу секущих
    static std::pair<var_type, numerical_result_code_t> next_argument_secant(
        const var_type& x1, const var_type& x2, const var_type& y1, const var_type& y2)
    {
        var_type x3 = (std::abs(y2) / (std::abs(y2) + std::abs(y1)) * x1 + std::abs(y1) / (std::abs(y2) + std::abs(y1)) * x2);
        // что-то пошло не так
        if (!std::isfinite(x3)) {
            return std::make_pair(std::numeric_limits<double>::quiet_NaN(), numerical_result_code_t::NumericalNanValues);
        }
        // попали в пределы
        if (x3 > x1 && x3 < x2) {
            return std::make_pair(x3, numerical_result_code_t::Converged);
        }
        return std::make_pair(std::numeric_limits<double>::quiet_NaN(), numerical_result_code_t::CustomCriteriaFailed);
    }
    /// @brief Расчет следующего приближения по методу бисекции
    static std::pair<var_type, numerical_result_code_t> next_argument_bisection(
        const var_type& x1, const var_type& x2, const var_type& y1, const var_type& y2)
    {
        var_type x3 = 0.5 * (x2 + x1);
        if (!std::isfinite(x3)) {
            return std::make_pair(std::numeric_limits<double>::quiet_NaN(), numerical_result_code_t::NumericalNanValues);
        }
        // тут мы не должны выпасть из интервала, можем попасть на границы ввиду точности и зациклится
        // reached_precision должно это исключить
        return std::make_pair(x3, numerical_result_code_t::Converged);
    }

    /// \brief Проверяет не пора ли завершать расчет (успешно или аварийно)
    ///  - сошлись по невзяке (успешно)
    ///  - получили NaN или Inf (аварийно)
    ///  В этих случаях также заполняет поля result_code и score в result.
    ///  В случае если analysis != nullptr :
    ///     -  в настройках указано сохранение истории аргументов,
    ///     то добавляет текущий аргумент в соответствующий массив analisis
    ///     -  в настройках указано сохранение истории невязока,
    ///     то добавляет текущую невязку в соответствующий массив analisis
    /// \param solver_parameters - параметры солвера
    /// \param r - текущее значение невязки
    /// \param argument - текущее значение аргумента
    /// \param analysis - указатель на служебную структуру, для хранения истории
    /// \param result - указатель на служебную структуру, для хранения резудьтата
    /// \return возвращает TRUE когда всё. ))
    static bool residual_exit_criterium(const fixed_bisectional_parameters_t& solver_parameters,
        var_type r, var_type argument, fixed_bisection_result_analysis_t<Dimension>* analysis,
        fixed_bisection_result_t<Dimension>* result)
    {
        if(solver_parameters.verbose){
            std::cerr<<"check:"<<r<<" with:"<<argument<<std::endl;
        }
        if(analysis != nullptr) {
            if(solver_parameters.argument_history)
                analysis->argument_history.push_back(argument);
            if(solver_parameters.residual_history)
                analysis->residual_history.push_back(r);
        }

        if(!std::isfinite(r)){
            result->result_code = numerical_result_code_t::NumericalNanValues;
            result->score = convergence_score_t::Error;
            return true;
        }

        if(std::fabs(r)<solver_parameters.residual_precision){
            result->result_code = numerical_result_code_t::Converged;
            result->score = convergence_score_t::Excellent;
            result->residuals = r;
            return true;
        }
        return false;
    }
    /// \brief поиск решения уравнения при полностью определённых входных параментах
    /// Критерии выхода 
    ///       1. По заданной абсолютной погрешности 
    ///       2. По количеству итераций на основе абсолютной погрешности. 
    ///       3. По машинному эпсилону
    /// Заметим, что 1 потенциально быстрее, чем 2 (т.к.Иллинойс), но может 
    /// не сработать при больших порядках чисел. Поэтому 2 нас подстраховывает.
    /// Проверка 3 нужна для больших аргументов
    /// \param residuals  Функция невязок
    /// \param result - указатель на служебную структуру, для хранения резудьтата
    /// \param x1 - нижняя граница аргумента
    /// \param x2 - верхняя граница аргумента
    /// \param x3 - первое приближение
    /// \param solver_parameters - параметры солвера
    /// \param analysis - указатель на служебную структуру, для хранения истории
    static void solve_limited(
            const fixed_bisectional_parameters_t& solver_parameters,
            fixed_system_t<Dimension>& residuals,
            var_type x1, var_type x2,
            fixed_bisection_result_t<Dimension>* result,
            fixed_bisection_result_analysis_t<Dimension>* analysis
            )
    {
        var_type& x3 = result->argument;

        var_type y1 = residuals.residuals(x1);
        var_type y2 = residuals.residuals(x2);
        var_type& y3 = result->residuals;
        y3 = residuals.residuals(x3);

        size_t& iteration = result->iteration_count;
        result->reached_precision = std::abs(x2 - x1);

        // Ниже в while проверка на диапазон abs(x2-x1) < argument_precision сделана через итерации, 
        // т.к. при малых значениях x1, x2 есть проблемы с такой проверкой
        result->max_allowed_iterations = get_max_allowed_iterations(
            std::abs(x2-x1),
            solver_parameters.argument_precision);


        int previous_residual_sign = 0; // знак функции на прошлой итерации
        for (iteration = 0; iteration <= result->max_allowed_iterations; ++iteration)
        {
            bool use_secant=false;
            std::tie(x3, result->result_code) = next_argument_value(solver_parameters, x1, x2, y1, y2, iteration,use_secant);
            if (result->result_code != numerical_result_code_t::Converged) {
                if (solver_parameters.verbose) {
                    std::cout << iteration << '\t' << result->max_allowed_iterations << std::endl;
                    std::cout << std::abs(x2 - x1) << '\t' << 2. * solver_parameters.argument_precision << '\t' << solver_parameters.residual_precision << std::endl;
                }
                result->score = convergence_score_t::Error;
                return;
            }
            y3 = residuals.residuals(x3);
            if (solver_parameters.verbose)
                std::cout <<std::setprecision(20)<< x3 << '\t' << y3 << '\t' << iteration << " of " << result->max_allowed_iterations << std::endl;
            if (residual_exit_criterium(solver_parameters, y3, x3, analysis, result)) {
                // Сошлись по невязке, она в приоритете
                result->score = convergence_score_t::Excellent;
                return;
            }

            if (y3 > 0) {
                x1 = x3;
                y1 = y3;
                if (use_secant && solver_parameters.use_Illinois && previous_residual_sign == +1) {
                    // Если в прошлый раз знак функции был тот же, можем сузить поиск
                    //https://en.wikipedia.org/wiki/Regula_falsi, раздел Illinois algorithm
                    y2 /= 2.0;
                }
                previous_residual_sign = +1;
            }
            else if (y3 < 0) {
                x2 = x3;
                y2 = y3;
                if (use_secant && solver_parameters.use_Illinois && previous_residual_sign == -1)
                    y1 /= 2.0;
                previous_residual_sign = -1;
            }
            // y3 == 0 не может быть, см. residual_exit_criterium

            // ЗНАЧЕНИЯ x1 или x2 изменидись, обновим достигнутую точность
            result->reached_precision = std::abs(x2 - x1);
            // диапазон поиска попал в абсолютная погрешность
            if (result->reached_precision <= solver_parameters.argument_precision)
                break;

            // диапазон поиска уже меньше машинного эпсилон x1 или x2)
            if (is_within_machine_epsilon(x1, x2))
                break;

        }

        // Если мы сюда попали, то не сошлись по невязке, она в приоритете
        // видимо задали недостижимую точнось

        // тут что-то пошло не так ))
        if (iteration >= result->max_allowed_iterations) {
            result->result_code = numerical_result_code_t::Converged;
            result->score = convergence_score_t::Satisfactory;
        }
        // иначе проанализируем, что у нас по точности решения
        else {
            result->result_code = numerical_result_code_t::Converged;
            if (result->reached_precision < solver_parameters.argument_precision) {
                // Достигнута точность по аргументу
                result->score = convergence_score_t::Excellent;
            }
            else {
                // Нельзя достигнуть точность по аргументу из-за точности double
                // Конкретнее - выскочили за машинный эпсилон, 
                // и тогда норма для абсолютной погрешности больше заданной в параметрах точности по аргументу.
                // Пример - если X был энтальпией, а хотели 0.00003 точность
                // ну что ж признаем результат "удовлетворительным"
                result->score = convergence_score_t::Good;
            }
        }
        return;
    }
public:
    /// @brief Запуск численного метода, начальное приближение рассчитывается на основе границ из solver_parameters
    /// @param residuals Функция невязок
    /// @param solver_parameters Настройки поиска
    /// @param result Результаты расчета
    /// @param analysis - указатель на служебную структуру, для хранения истории
    static void solve(
        const fixed_bisectional_parameters_t& solver_parameters,
        fixed_system_t<Dimension>& residuals,
        fixed_bisection_result_t<Dimension>* result,
        fixed_bisection_result_analysis_t<Dimension>* analysis = nullptr)
    {
        var_type initial_argument = std::numeric_limits<double>::quiet_NaN();
        solve(solver_parameters, initial_argument, residuals, result, analysis);
    }
    /// @brief Запуск численного метода
    /// @param residuals Функция невязок
    /// @param initial_argument Начальное приближение
    /// @param solver_parameters Настройки поиска
    /// @param result Результаты расчета
    /// @param analysis - указатель на служебную структуру, для хранения истории
    static void solve(
        const fixed_bisectional_parameters_t& solver_parameters,
        var_type initial_argument,
        fixed_system_t<Dimension>& residuals,
        fixed_bisection_result_t<Dimension>* result,
        fixed_bisection_result_analysis_t<Dimension>* analysis = nullptr
    )
    {
        // Проверка корректности границ по x
        var_type minx = solver_parameters.argument_limit_min;
        var_type maxx = solver_parameters.argument_limit_max;
        if(!std::isfinite(minx)){
            throw std::logic_error("Не задана нижняя граница");
        }
        if(!std::isfinite(maxx)){
            throw std::logic_error("Не задана верхняя граница");
        }

        // Подготовка и проверка начального приближения 
        function_type& argument = result->argument;
        if (!std::isfinite(initial_argument)) {
            argument = (maxx + minx) / 2.0;
        }
        else {
            if (argument < minx || argument > maxx) {
                throw std::logic_error("Не верно задано начальное приближение");
            }
            argument = initial_argument;
        }

        // Проверка сходимости по невязкам на границах интервала поиска
        // todo: (подумал - это неверно) Раз тут уже считаются невязки, значения miny, minx стоит передавать в solve_limited
        if(solver_parameters.check_boundary_before){
            var_type miny = residuals.residuals(minx);
            var_type maxy = residuals.residuals(maxx);
            if (residual_exit_criterium(solver_parameters, miny, minx, analysis, result))
                return;
            if (residual_exit_criterium(solver_parameters, maxy, maxx, analysis, result))
                return;
        }

        // Проверка сходимости по невязкам в точке начального приближения
        var_type& r = result->residuals;
        r = residuals.residuals(argument);
        if (residual_exit_criterium(solver_parameters, r, argument, analysis, result)) 
            return;

        if (r > 0)
            minx = argument;
        if (r < 0)
            maxx = argument;

        // Запуск расчета методом бисекции
        result->result_code = numerical_result_code_t::NotConverged;
        result->score = convergence_score_t::Excellent;
        
        // Здесь подготовленная переменная argument передается в составе result
        solve_limited(solver_parameters, residuals, minx, maxx,
            result, analysis);
    }

};
