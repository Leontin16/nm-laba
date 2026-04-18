#include "algorithm.h"


// построение равномерной сетки из n+1 узлов на отрезке [a, b]
std::vector<double> buildUniformGrid(double a, double b, int n) {
    std::vector<double> x(n + 1);
    double h = (b - a) / n;
    for (int i = 0; i <= n; i++) {
        x[i] = a + i * h;
    }
    return x; 
}


std::vector<double> computeSteps(const std::vector<double>& x) {
    int n = (int)x.size() - 1;
    std::vector<double> h(n);
    for (int i = 0; i < n; i++) {
        h[i] = x[i + 1] - x[i];
    }
    return h;
}


// формирование 3-диагональной СЛАУ для метода прогонки
TridiagonalSystem buildSplineSystem(const std::vector<double>& x, 
                                    const std::vector<double>& f, 
                                    double bc_left, 
                                    double bc_right) {
    int n = (int)x.size() - 1; // число участков
    int m = n - 1;             // число внутренних узлов

    std::vector<double> h = computeSteps(x);

    TridiagonalSystem TDS;
    TDS.lower.resize(m, 0.0);
    TDS.main.resize(m, 0.0);
    TDS.upper.resize(m, 0.0);
    TDS.rhs.resize(m, 0.0);

    for (int i = 0; i < m; i++) {
        int k = i + 1; // индекс внутреннего узла

        // коэффициент при c[k-1]
        TDS.lower[i] = h[k - 1];

        // коэффициент при c[k]
        TDS.main[i] = 2.0 * (h[k - 1] + h[k]);

        // коэффициент при c[k+1]
        TDS.upper[i] = h[k];

        // правая часть
        TDS.rhs[i] = 6.0 * ((f[k + 1] - f[k]) / h[k] - (f[k] - f[k - 1]) / h[k - 1]);

        if (i == 0) {
            // левое граничное условие 
            TDS.rhs[i] -= h[k - 1] * bc_left;
        }

        if (i == m - 1) {
            TDS.rhs[i] -= h[k] * bc_right;
        }
    }

    TDS.lower[0] = 0.0;
    TDS.upper[m - 1] = 0.0;

    return TDS;
}