#include <iostream>
#include <cmath>
#include <vector>
#include "algorithm.h"

// --- Функции из твоего изначального файла (для варианта a) ---
double phi(double x) {
    if (x <= 0) return x * x * x + 3.0 * x * x;
    return -x * x * x + 3.0 * x * x;
}

double phi_double_derivative(double x) {
    if (x <= 0) return 6.0 * x + 6.0;
    return -6.0 * x + 6.0;
}

// =====================================================================
// ТЕСТ 1: Проверка построения сетки и вычисления шагов
// =====================================================================
bool test_grid_and_steps() {
    std::cout << "[TEST 1] Вектор сетки и шагов... ";
    auto x = buildUniformGrid(0.0, 1.0, 2);
    if (x.size() != 3 || std::abs(x[0] - 0.0) > 1e-9 || std::abs(x[1] - 0.5) > 1e-9 || std::abs(x[2] - 1.0) > 1e-9) {
        std::cout << "FAIL (ошибка сетки)\n";
        return false;
    }

    auto h = computeSteps(x);
    if (h.size() != 2 || std::abs(h[0] - 0.5) > 1e-9 || std::abs(h[1] - 0.5) > 1e-9) {
        std::cout << "FAIL (ошибка шагов)\n";
        return false;
    }
    std::cout << "OK\n";
    return true;
}

// =====================================================================
// ТЕСТ 2: Проверка решателя 3-диагональной СЛАУ (Метод прогонки)
// =====================================================================
bool test_tridiagonal_solver() {
    std::cout << "[TEST 2] Решение 3-диагональной СЛАУ... ";
    TridiagonalSystem sys;
    sys.lower = {0.0, -1.0, -1.0}; // Первый элемент не используется
    sys.main  = {2.0,  2.0,  2.0};
    sys.upper = {-1.0, -1.0, 0.0}; // Последний элемент не используется
    sys.rhs   = {1.0,  0.0,  1.0};
    
    // Решение этой системы должно быть {1.0, 1.0, 1.0}
    auto res = solveTridiagonalSystem(sys);
    for (double val : res) {
        if (std::abs(val - 1.0) > 1e-9) {
            std::cout << "FAIL (ошибка прогонки)\n";
            return false;
        }
    }
    std::cout << "OK\n";
    return true;
}

// =====================================================================
// ТЕСТ 3: Проверка формирования системы (твой изначальный тест)
// =====================================================================
bool test_spline_system_build() {
    std::cout << "[TEST 3] Сборка СЛАУ для сплайна (phi)... \n";
    double a = -1.0, b = 1.0;
    int n = 4;

    std::vector<double> x = buildUniformGrid(a, b, n);
    std::vector<double> f(n + 1);
    for (int i = 0; i <= n; i++) f[i] = phi(x[i]);

    double bc_left = phi_double_derivative(a);
    double bc_right = phi_double_derivative(b);

    TridiagonalSystem TDS = buildSplineSystem(x, f, bc_left, bc_right);
    int m = n - 1; 

    bool ok = true;
    for (int i = 0; i < m; i++) {
        int k = i + 1;
        double c_prev = (i == 0) ? bc_left : phi_double_derivative(x[k - 1]);
        double c_curr = phi_double_derivative(x[k]);
        double c_next = (i == m - 1) ? bc_right : phi_double_derivative(x[k + 1]);

        double lhs = TDS.lower[i] * c_prev + TDS.main[i] * c_curr + TDS.upper[i] * c_next;
        double diff = std::abs(lhs - TDS.rhs[i]);

        if (diff > 1e-10) ok = false;
    }
    if (ok) std::cout << "         -> OK\n";
    else std::cout << "         -> FAIL (невязка превышает допуск)\n";
    
    return ok;
}

// =====================================================================
// ТЕСТ 4: Проверка коэффициентов и интерполяции в узлах
// =====================================================================
bool test_spline_interpolation_at_nodes() {
    std::cout << "[TEST 4] Интерполяция в узловых точках (Вариант 2)... ";
    int n = 5;
    auto x = buildUniformGrid(0.0, 1.0, n);
    std::vector<double> f(n + 1);
    for (int i = 0; i <= n; i++) f[i] = main_func(x[i]);

    // Естественные граничные условия
    auto segments = computeSplineCoefficients(x, f, 0.0, 0.0);

    for (int i = 0; i <= n; i++) {
        double s_val = evaluateSpline(segments, x[i]);
        if (std::abs(s_val - f[i]) > 1e-9) {
            std::cout << "FAIL (Сплайн не проходит через узел x=" << x[i] << ")\n";
            return false;
        }
    }
    std::cout << "OK\n";
    return true;
}

// =====================================================================
// ТЕСТ 5: Проверка функции вычисления норм погрешностей
// =====================================================================
bool test_error_norms() {
    std::cout << "[TEST 5] Вычисление норм погрешностей... ";
    int n = 4;
    auto x = buildUniformGrid(0.0, 1.0, n);
    std::vector<double> f(n + 1);
    for (int i = 0; i <= n; i++) f[i] = main_func(x[i]);

    auto segments = computeSplineCoefficients(x, f, 0.0, 0.0);
    
    // Проверяем на густой сетке (10 * n)
    ErrorNorms norms = calculateErrorNorms(segments, main_func, main_func_d1, main_func_d2, 0.0, 1.0, n * 10);
    
    // Погрешность должна быть адекватной (не NaN, не бесконечность и больше нуля)
    if (norms.max_f < 0.0 || std::isnan(norms.max_f) || norms.max_f > 1.0) {
        std::cout << "FAIL (Некорректная норма погрешности)\n";
        return false;
    }
    std::cout << "OK (max_f = " << norms.max_f << ")\n";
    return true;
}

// =====================================================================
// Главная точка входа
// =====================================================================
int main() {
    std::cout << "=== RUNNING TESTS (Total: 5) ===\n\n";

    int passed = 0;
    if (test_grid_and_steps()) { std::cout << "[TEST 1] Grid and Steps... OK\n"; passed++; }
    if (test_tridiagonal_solver()) { std::cout << "[TEST 2] Tridiagonal Solver... OK\n"; passed++; }
    
    // Внутри test_spline_system_build я уже добавил вывод, так что просто вызываем
    if (test_spline_system_build()) passed++;

    if (test_spline_interpolation_at_nodes()) { std::cout << "[TEST 4] Interpolation at nodes... OK\n"; passed++; }
    if (test_error_norms()) { std::cout << "[TEST 5] Error norms calculation... OK\n"; passed++; }

    std::cout << "\n=== RESULT: " << passed << " / 5 TESTS PASSED ===\n";
    
    return (passed == 5) ? 0 : 1;
}