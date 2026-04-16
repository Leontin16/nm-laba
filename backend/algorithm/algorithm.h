#pragma once
#include <vector>


struct TridiagonalSystem {
    std::vector<double> lower; // нижнаяя диагональ
    std::vector<double> main; // главная диагональ
    std::vector<double> upper; // верхняя диагональ
    std::vector<double> rhs; // правая часть
};


// построение равномерной сетки из n+1 узлов на отрезке [a, b]
std::vector<double> buildUniformGrid(double a, double b, int n);


// вычисление шага сетки
std::vector<double> computeSteps(const std::vector<double>& x);


// формирование 3-диагональной СЛАУ для метода прогонки
TridiagonalSystem buildSplineSystem(const std::vector<double>& x, 
                                    const std::vector<double>& f, 
                                    double bc_left = 0.0, 
                                    double bc_right = 0.0);
