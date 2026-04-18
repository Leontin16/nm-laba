#include <iostream>
#include <cmath>
#include <vector>
#include "algorithm.h"


// Тестовая функция φ(x) из задания (вариант a)
// φ(x) = x^3 + 3x^2,  x in [-1, 0]
// φ(x) = -x^3 + 3x^2, x in [0, 1] 
double phi(double x) {
    if (x <= 0) {
        return x*x*x + 3.0*x*x;
    } else {
        return -x*x*x + 3.0*x*x;
    }
}


// Точная вторая производная φ''(x)
// φ''(x) = 6x + 6,  x in [-1, 0]
// φ''(x) = -6x + 6, x in [0, 1]
double phi_double_derivative(double x) {
    if (x <= 0) {
        return 6.0*x + 6.0;
    } else {
        return -6.0*x + 6.0;
    }
}


int main() {
    // параметры сетки
    double a = -1.0;
    double b = 1.0;
    int n = 4;

    // построение сетки и вычисление значений второй производной в узлах
    std::vector<double> x = buildUniformGrid(a, b, n);
    std::vector<double> f(n + 1);
    for (int i = 0; i <= n; i++) {
        f[i] = phi(x[i]);
    }

    // граничные условия: S''(a) = φ''(-1), S''(b) = φ''(1)
    double bc_left = phi_double_derivative(a); // φ''(-1) = 0
    double bc_right = phi_double_derivative(b); // φ''(1) = 0

    std::cout << "bc_left = " << bc_left << "\n";
    std::cout << "bc_right = " << bc_right << "\n";

    // построение СЛАУ
    TridiagonalSystem TDS = buildSplineSystem(x, f, bc_left, bc_right);

    int m = n - 1; 

    std::cout << "Tridiagonal system (lower | main | upper | rhs):\n";
    for (int i = 0; i < m; i++) {
        std::cout << " [" << i << "] " <<TDS.lower[i] << " | " << TDS.main[i] << " | " << TDS.upper[i] << " | " << TDS.rhs[i] << "\n";
    }

    bool ok = true;
    double h = (b - a) / n;
    for (int i = 0; i < m; i++) {
        int k = i + 1;
        double c_prev = (i == 0) ? bc_left : phi_double_derivative(x[k - 1]);
        double c_curr = phi_double_derivative(x[k]);
        double c_next = (i == m - 1) ? bc_right : phi_double_derivative(x[k + 1]);

        double lhs = TDS.lower[i] * c_prev + TDS.main[i] * c_curr + TDS.upper[i] * c_next;
        double diff = std::abs(lhs - TDS.rhs[i]);

        std::cout << " equation " << i << ": lhs = " << lhs << ", rhs = " << TDS.rhs[i] << ", |diff| = " << diff << "\n";

        if (diff < 1e-10) {
            std::cout << "  OK\n";
        } else {
            std::cout << "  FAIL\n";
            ok = false;
        }
    }

    std::cout << "\nTest " << (ok ? "PASSED" : "FAILED") << "\n";
    return ok ? 0 : 1;
}