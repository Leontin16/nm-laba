// accuracy_engine_simple.cpp
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include "algorithm.h"

// ------------------------------------------------------------
// Вспомогательные функции
// ------------------------------------------------------------
std::vector<double> buildControlGrid(double a, double b, int N) {
    std::vector<double> x(N + 1);
    double h = (b - a) / N;
    for (int i = 0; i <= N; ++i)
        x[i] = a + i * h;
    return x;
}

struct ErrorSummary {
    double max_f, x_max_f;
    double max_d1, x_max_d1;
    double max_d2, x_max_d2;
};

void writeTable1(std::ofstream& out, const std::vector<SplineSegment>& seg) {
    out << "i,x_{i-1},x_i,a,b,c,d\n";
    for (size_t i = 0; i < seg.size(); ++i) {
        const auto& s = seg[i];
        out << (i+1) << "," << s.x_prev << "," << s.x_curr << ","
            << s.a << "," << s.b << "," << s.c << "," << s.d << "\n";
    }
}

void writeTable2(std::ofstream& out,
                 const std::vector<SplineSegment>& seg,
                 const std::vector<double>& x_control,
                 double (*ref_f)(double),
                 double (*ref_d1)(double),
                 double (*ref_d2)(double),
                 ErrorSummary& sum)
{
    out << "j,x,F(x),S(x),diff_F,F'(x),S'(x),diff_F',F''(x),S''(x),diff_F''\n";
    sum = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    for (size_t j = 0; j < x_control.size(); ++j) {
        double x = x_control[j];
        double Fx   = ref_f(x);
        double Sx   = evaluateSpline(seg, x);
        double dF   = Fx - Sx;
        double Fd1  = ref_d1(x);
        double Sd1  = evaluateSplineDeriv1(seg, x);
        double dD1  = Fd1 - Sd1;
        double Fd2  = ref_d2(x);
        double Sd2  = evaluateSplineDeriv2(seg, x);
        double dD2  = Fd2 - Sd2;

        out << j << "," << x << ","
            << Fx << "," << Sx << "," << dF << ","
            << Fd1 << "," << Sd1 << "," << dD1 << ","
            << Fd2 << "," << Sd2 << "," << dD2 << "\n";

        double af = std::abs(dF);
        if (af > sum.max_f) { sum.max_f = af; sum.x_max_f = x; }
        double ad1 = std::abs(dD1);
        if (ad1 > sum.max_d1) { sum.max_d1 = ad1; sum.x_max_d1 = x; }
        double ad2 = std::abs(dD2);
        if (ad2 > sum.max_d2) { sum.max_d2 = ad2; sum.x_max_d2 = x; }
    }
}

void writeErrorSummary(std::ofstream& out, const ErrorSummary& sum) {
    out << "=== ERROR SUMMARY ===\n";
    out << "max |F - S|      = " << sum.max_f  << " at x = " << sum.x_max_f  << "\n";
    out << "max |F'- S'|     = " << sum.max_d1 << " at x = " << sum.x_max_d1 << "\n";
    out << "max |F''- S''|   = " << sum.max_d2 << " at x = " << sum.x_max_d2 << "\n";
}

void processFunction(const std::string& name,
                     double a, double b,
                     double bc_left, double bc_right,
                     double (*ref_f)(double),
                     double (*ref_d1)(double),
                     double (*ref_d2)(double),
                     int n, int N_control)
{
    // Узлы основной сетки
    auto x_nodes = buildUniformGrid(a, b, n);
    std::vector<double> f_nodes(x_nodes.size());
    for (size_t i = 0; i < x_nodes.size(); ++i)
        f_nodes[i] = ref_f(x_nodes[i]);

    // Построение сплайна
    auto segments = computeSplineCoefficients(x_nodes, f_nodes, bc_left, bc_right);

    // Контрольная сетка
    auto x_control = buildControlGrid(a, b, N_control);

    // Сохраняем отчёты в файлы
    std::string prefix = name + "_n" + std::to_string(n);
    std::ofstream f1(prefix + "_table1.csv");
    writeTable1(f1, segments);
    f1.close();

    std::ofstream f2(prefix + "_table2.csv");
    ErrorSummary sum;
    writeTable2(f2, segments, x_control, ref_f, ref_d1, ref_d2, sum);
    f2.close();

    std::ofstream fsum(prefix + "_errors.txt");
    writeErrorSummary(fsum, sum);
    fsum.close();

    std::cout << "Saved: " << prefix << "_*.csv/txt\n";
}

void convergenceStudy(const std::string& name,
                      double a, double b,
                      double bc_left, double bc_right,
                      double (*ref_f)(double),
                      double (*ref_d1)(double),
                      double (*ref_d2)(double),
                      const std::vector<int>& n_vals,
                      int control_mult = 10)
{
    std::string fname = name + "_convergence.csv";
    std::ofstream out(fname);
    out << "n,max_err_F,max_err_F',max_err_F''\n";

    std::vector<double> prev_n, prev_err_f, prev_err_d1, prev_err_d2;
    for (int n : n_vals) {
        auto x_nodes = buildUniformGrid(a, b, n);
        std::vector<double> f_nodes(x_nodes.size());
        for (size_t i = 0; i < x_nodes.size(); ++i) f_nodes[i] = ref_f(x_nodes[i]);

        auto seg = computeSplineCoefficients(x_nodes, f_nodes, bc_left, bc_right);
        auto x_control = buildControlGrid(a, b, n * control_mult);

        double max_f = 0, max_d1 = 0, max_d2 = 0;
        for (double x : x_control) {
            double ef = std::abs(ref_f(x) - evaluateSpline(seg, x));
            double ed1 = std::abs(ref_d1(x) - evaluateSplineDeriv1(seg, x));
            double ed2 = std::abs(ref_d2(x) - evaluateSplineDeriv2(seg, x));
            if (ef > max_f) max_f = ef;
            if (ed1 > max_d1) max_d1 = ed1;
            if (ed2 > max_d2) max_d2 = ed2;
        }
        out << n << "," << max_f << "," << max_d1 << "," << max_d2 << "\n";
    }
    out.close();
    std::cout << "Convergence data saved to " << fname << "\n";
}

// ------------------------------------------------------------
// main
// ------------------------------------------------------------
int main() {
    // Параметры варианта 2
    double a = 0.0, b = 1.0;
    double bc_left = 0.0, bc_right = 0.0;

    // Часть b) основная функция f(x)
    processFunction("f", a, b, bc_left, bc_right,
                    main_func, main_func_d1, main_func_d2,
                    20, 200);

    convergenceStudy("f", a, b, bc_left, bc_right,
                     main_func, main_func_d1, main_func_d2,
                     {5, 10, 20, 40, 80, 160});

    // Часть c) осциллирующая функция
    processFunction("osc", a, b, bc_left, bc_right,
                    osc_func, osc_func_d1, osc_func_d2,
                    20, 200);

    convergenceStudy("osc", a, b, bc_left, bc_right,
                     osc_func, osc_func_d1, osc_func_d2,
                     {5, 10, 20, 40, 80, 160});

    // Тестовая функция phi(x)
    processFunction("phi", -1.0, 1.0,
                    test_phi_d2(-1.0), test_phi_d2(1.0),
                    test_phi, test_phi_d1, test_phi_d2,
                    4, 40);

    std::cout << "All calculations finished.\n";
    return 0;
}