#include "algorithm.h"
#include <stdexcept>

std::vector<double> solveTridiagonalSystem(const TridiagonalSystem& sys) {
    int m = (int)sys.rhs.size();
    if (m == 0) return {};

    std::vector<double> alpha(m);
    std::vector<double> beta(m);
    std::vector<double> x_res(m);

    // Прямой ход прогонки
    // Для первого уравнения: b0*x0 + c0*x1 = f0
    alpha[0] = -sys.upper[0] / sys.main[0];
    beta[0] = sys.rhs[0] / sys.main[0];

    for (int i = 1; i < m; i++) {
        double denominator = sys.main[i] + sys.lower[i] * alpha[i - 1];
        if (i < m - 1) {
            alpha[i] = -sys.upper[i] / denominator;
        }
        beta[i] = (sys.rhs[i] - sys.lower[i] * beta[i - 1]) / denominator;
    }

    // Обратный ход
    x_res[m - 1] = beta[m - 1];
    for (int i = m - 2; i >= 0; i--) {
        x_res[i] = alpha[i] * x_res[i + 1] + beta[i];
    }

    return x_res;
}

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

TridiagonalSystem buildSplineSystem(const std::vector<double>& x, 
                                    const std::vector<double>& f, 
                                    double bc_left, 
                                    double bc_right) {
    int n = (int)x.size() - 1; 
    int m = n - 1; // Число внутренних узлов, где ищутся c_i [cite: 29]

    std::vector<double> h = computeSteps(x);
    TridiagonalSystem TDS;
    TDS.lower.resize(m);
    TDS.main.resize(m);
    TDS.upper.resize(m);
    TDS.rhs.resize(m);

    for (int i = 0; i < m; i++) {
        int k = i + 1; // Узел сетки x[1]...x[n-1]

        TDS.lower[i] = h[k - 1];
        TDS.main[i] = 2.0 * (h[k - 1] + h[k]);
        TDS.upper[i] = h[k];
        
        // Правая часть на основе разделенных разностей [cite: 3]
        TDS.rhs[i] = 6.0 * ((f[k + 1] - f[k]) / h[k] - (f[k] - f[k - 1]) / h[k - 1]);

        // Естественные граничные условия S''(a)=bc_left, S''(b)=bc_right [cite: 13, 14]
        if (i == 0) TDS.rhs[i] -= h[k - 1] * bc_left;
        if (i == m - 1) TDS.rhs[i] -= h[k] * bc_right;
    }

    if (m > 0) {
        TDS.lower[0] = 0.0;
        TDS.upper[m - 1] = 0.0;
    }

    return TDS;
}

std::vector<SplineSegment> computeSplineCoefficients(
    const std::vector<double>& x, 
    const std::vector<double>& f, 
    double bc_left, 
    double bc_right) 
{
    int n = (int)x.size() - 1;
    std::vector<double> h = computeSteps(x);

    // 1. Формируем и решаем систему для внутренних узлов c_1 ... c_{n-1}
    TridiagonalSystem sys = buildSplineSystem(x, f, bc_left, bc_right);
    std::vector<double> c_internal = solveTridiagonalSystem(sys);

    // 2. Формируем полный вектор c (от c_0 до c_n)
    // По заданию S''(a) = bc_left и S''(b) = bc_right [cite: 13]
    std::vector<double> c_full(n + 1);
    c_full[0] = bc_left;
    c_full[n] = bc_right;
    for (int i = 0; i < (int)c_internal.size(); ++i) {
        c_full[i + 1] = c_internal[i];
    }

    // 3. Вычисляем a, b, d для каждого участка i = 1...n
    std::vector<SplineSegment> segments(n);
    for (int i = 1; i <= n; ++i) {
        int idx = i - 1; // индекс в массиве h и segments
        
        segments[idx].x_prev = x[i - 1];
        segments[idx].x_curr = x[i];
        
        // Согласно канонической записи S(x_i) = a_i 
        segments[idx].a = f[i];
        segments[idx].c = c_full[i];
        segments[idx].d = (c_full[i] - c_full[i - 1]) / h[idx];
        segments[idx].b = (f[i] - f[i - 1]) / h[idx] + (2.0 * c_full[i] + c_full[i - 1]) * h[idx] / 6.0;
    }

    return segments;
}


double evaluateSpline(const std::vector<SplineSegment>& segments, double x) {
    for (const auto& seg : segments) {
        if (x >= seg.x_prev && x <= seg.x_curr) {
            double dx = x - seg.x_curr;
            return seg.a + seg.b * dx + (seg.c / 2.0) * dx * dx + (seg.d / 6.0) * dx * dx * dx;
        }
    }
    return 0.0; // Точка вне диапазона
}

double evaluateSplineDeriv1(const std::vector<SplineSegment>& segments, double x) {
    for (const auto& seg : segments) {
        if (x >= seg.x_prev && x <= seg.x_curr) {
            double dx = x - seg.x_curr;
            return seg.b + seg.c * dx + (seg.d / 2.0) * dx * dx;
        }
    }
    return 0.0;
}

double evaluateSplineDeriv2(const std::vector<SplineSegment>& segments, double x) {
    for (const auto& seg : segments) {
        if (x >= seg.x_prev && x <= seg.x_curr) {
            double dx = x - seg.x_curr;
            return seg.c + seg.d * dx;
        }
    }
    return 0.0;
}

#include <cmath>

double test_phi(double x) {
    if (x <= 0) return std::pow(x, 3) + 3 * std::pow(x, 2);
    return -std::pow(x, 3) + 3 * std::pow(x, 2);
}

double test_phi_d1(double x) {
    if (x <= 0) return 3 * std::pow(x, 2) + 6 * x;
    return -3 * std::pow(x, 2) + 6 * x;
}

double test_phi_d2(double x) {
    if (x <= 0) return 6 * x + 6;
    return -6 * x + 6;
}

double main_func(double x) {
    return std::pow(1.0 + x * x, 1.0 / 3.0);
}

// Производные для варианта 2 (вычислены аналитически)
double main_func_d1(double x) {
    return (2.0 * x) / (3.0 * std::pow(1.0 + x * x, 2.0 / 3.0));
}

double main_func_d2(double x) {
    double val = 1.0 + x * x;
    return (2.0 / (3.0 * std::pow(val, 2.0 / 3.0))) - (8.0 * x * x / (9.0 * std::pow(val, 5.0 / 3.0)));
}

// Добавь реализацию осциллирующей функции (пункт 'c' задания)
double osc_func(double x) {
    return main_func(x) + std::cos(10.0 * x);
}

// Производные для осциллирующей функции
double osc_func_d1(double x) {
    return main_func_d1(x) - 10.0 * std::sin(10.0 * x);
}

double osc_func_d2(double x) {
    return main_func_d2(x) - 100.0 * std::cos(10.0 * x);
}

// Функция для расчета погрешностей по Таблице 5
ErrorNorms calculateErrorNorms(const std::vector<SplineSegment>& segments, 
                               double (*ref_f)(double), 
                               double (*ref_d1)(double), 
                               double (*ref_d2)(double),
                               double a, double b, int n_control) 
{
    ErrorNorms norms = {0.0, 0.0, 0.0};
    double h_control = (b - a) / n_control;

    for (int i = 0; i <= n_control; ++i) {
        double x = a + i * h_control;
        
        double err_f = std::abs(ref_f(x) - evaluateSpline(segments, x));
        double err_d1 = std::abs(ref_d1(x) - evaluateSplineDeriv1(segments, x));
        double err_d2 = std::abs(ref_d2(x) - evaluateSplineDeriv2(segments, x));

        if (err_f > norms.max_f) norms.max_f = err_f;
        if (err_d1 > norms.max_d1) norms.max_d1 = err_d1;
        if (err_d2 > norms.max_d2) norms.max_d2 = err_d2;
    }
    return norms;
}