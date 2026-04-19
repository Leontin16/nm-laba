#pragma once
#include <vector>

// Структура для хранения коэффициентов 3-диагональной системы
struct TridiagonalSystem {
    std::vector<double> lower; // нижняя диагональ (a_i)
    std::vector<double> main;  // главная диагональ (b_i)
    std::vector<double> upper; // верхняя диагональ (c_i)
    std::vector<double> rhs;   // правая часть (f_i)
};

// Решение 3-диагональной СЛАУ методом прогонки
std::vector<double> solveTridiagonalSystem(const TridiagonalSystem& sys);

// Построение равномерной сетки
std::vector<double> buildUniformGrid(double a, double b, int n);

// Вычисление шага сетки
std::vector<double> computeSteps(const std::vector<double>& x);

// Структура для одного сегмента сплайна на интервале [x_{i-1}, x_i]
struct SplineSegment {
    double x_prev; // x_{i-1}
    double x_curr; // x_i
    double a, b, c, d; // Коэффициенты
};

// Итоговая функция: вычисляет все коэффициенты для n участков
std::vector<SplineSegment> computeSplineCoefficients(
    const std::vector<double>& x, 
    const std::vector<double>& f, 
    double bc_left = 0.0, 
    double bc_right = 0.0
);

// Формирование СЛАУ для сплайна
TridiagonalSystem buildSplineSystem(const std::vector<double>& x, 
                                    const std::vector<double>& f, 
                                    double bc_left = 0.0, 
                                    double bc_right = 0.0); 

// Вычисление значения сплайна и его производных в точке x
double evaluateSpline(const std::vector<SplineSegment>& segments, double x);
double evaluateSplineDeriv1(const std::vector<SplineSegment>& segments, double x);
double evaluateSplineDeriv2(const std::vector<SplineSegment>& segments, double x);

// Структура для хранения погрешностей (для Таблицы 5 из задания)
struct ErrorNorms {
    double max_f;  // по функции
    double max_d1; // по 1-й производной
    double max_d2; // по 2-й производной
};                                    

// Тестовая функция и её производные
double test_phi(double x);
double test_phi_d1(double x);
double test_phi_d2(double x);

// Основная функция (Вариант 2: f(x) = (1+x^2)^(1/3))
double main_func(double x);
double main_func_d1(double x);
double main_func_d2(double x);


double osc_func(double x);
double osc_func_d1(double x);
double osc_func_d2(double x);

ErrorNorms calculateErrorNorms(const std::vector<SplineSegment>& segments, 
                               double (*ref_f)(double), 
                               double (*ref_d1)(double), 
                               double (*ref_d2)(double),
                               double a, double b, int n_control);