# Решатель Ньютона - Рафсона для решения систем нелинейных уравнений фиксированный размерности

Данный солвер решает систему нелинейных уравнений заранее известной (фиксированной) размерности с помощью алгоритма Ньютона-Рафсона.
Особенностью данного солвера является учет фиксированности размерности, который позволяет не использовать динамическую память для хранения векторов.

## Оглавление
* [Назначение библиотеки](#Назначение-библиотеки)
  * [Метод Ньютона - Рафсона](#Метод-Ньютона---Рафсона)
  * [Пример ипользования](#Пример-использования)
* [Инициализация системы нелинейных уравнений (fixed_system.h)](#fixed_system.h)
* [Решатель Ньютона - Рафсона (fixed_nonlinear_solver.h)](#fixed_nonlinear_solver.h)
* [Пример использования библиотеки](#Пример-использования-библиотеки)


## Назначение библиотеки

### Метод Ньютона - Рафсона

Метод Ньютона - Рафсона - это итерационный чисенный метод нахождения корня уравнения (системы уравнений). Поиск решения осуществляется путём построения последовательных приближений и основан на принципах простой итерации.

### Пример использования

Рассмотрим следующую систему уравнений: 

```math
f(n) =
  \begin{cases}
    n/2       & \quad \text{if } n \text{ is even}\\
    -(n+1)/2  & \quad \text{if } n \text{ is odd}
  \end{cases}
```


<div align="center">
</center><h2>Структура библиотеки</h2></center>
</div>

<pre>
├───fixed
│   │   array_ext.h - Дополнительные матричные операции
│   │   fixed.h - Этот файл надо подключить для использований функций библиотеки
│   │   fixed_linear_solver.h - Решатель системы линейных уравнений
│   │   fixed_nonlinear_solver.h - Решатель системы НЕлинейных уравнений
│   │   fixed_system.h - Специализации систем различной размерности
│   │
│   └───line_search - Алгоритмы линейного поиска
│           divider.h - Линейный поиск методом дробления шага
│           golden_section.h - Линейный поиск методом золотого сечения
</pre>

<div align="center">
</center><h2>fixed_system.h</h2></center>
</div>

Данный *.h* файл представляет собой структуры с помощью которых можно задать системы нелинейных уравнений фиксированной размерности с помощью массиовов, что позволяет не использовать динамическую память для хранения векторов.

<div align="center">
  
|Тип системы нелинейных уравнений|Размерность|
|:----:|:----------:|
|Скалярный случай|<1>|
|Векторный случай|<*n*>|

</div>

Также здесь реализованы:
1. Численный расчет Якобиана для скалярного и векторного случая 
2. Расчет целевой функции методом суммы квадратов 

<div align="center">
</center><h2>fixed_nonlinear_solver.h</h2></center>
</div>

В данном *.h* файле реализован алгоритм Ньютона - Рафсона. С помощью которого решается заданная система нелинейных уравнений. Код представляет выбор метода линейного поиска:
1. Линейный поиск методом золотого сечения 
2. Линейный поиск методом дробления шага 

Настройка, анализ и результат работы алгоритма реализован с помощью следующих структур:

<div align="center">
  
|Название структуры|Свойства|
|:----:|:----------:|
|<fixed_solver_parameters_t>|В данной структуре задаются настроечные параметры алгоритма|
|<fixed_solver_result_t>|В данную структуру записывается результат работы решателя|
|<fixed_solver_result_analysis_t>|В данную структуру записывается анализ работы решателя|
|<fixed_solver_constraints>|В данной структуре задаются ограничения для решателя|

</div>
  
<div align="center">
</center><h2>Пример использования библиотеки</h2></center>
</div>
  
```C++
TEST(NewtonRaphson, UseCase)
{
    // Класс, для системы размерности <2> - Векторный случай
    class sample_system : public fixed_system_t<2>
    {
        using fixed_system_t<2>::var_type;
    public:
        // Задание функции невязок
        var_type residuals(const var_type& x) {
            return
            {
                pow(x[0] - 2, 3),
                pow(x[1] - 1, 3)
            };
        }
    };
   
    sample_system test;
    fixed_solver_parameters_t<2, 0> parameters;
    fixed_solver_result_t<2> result;
    // Решение системы нелинейныйх уравнений <2> с помощью решателя Ньютона - Рафсона
    fixed_newton_raphson<2>::solve_dense(test, { 0, 0 }, parameters, & result);

}
```
