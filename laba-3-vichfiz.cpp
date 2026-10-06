#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <utility>
#include "gnuplot-iostream.h"

double f(double x) {
    return 1.0 / (1.0 + 25.0 * x * x);
}

double df_exact(double x) {
    return -50.0 * x / std::pow(1.0 + 25.0 * x * x, 2);
}

double d2f_exact(double x) {
    double denom = std::pow(1.0 + 25.0 * x * x, 3);
    return (3750.0 * x * x - 50.0) / denom;
}


double calc_diff1_m1_n3(double x1, double h) {
    double y0 = f(x1 - h);
    double y1 = f(x1);
    double y2 = f(x1 + h);
    double y3 = f(x1 + 2.0 * h);
    return (-2.0 * y0 - 3.0 * y1 + 6.0 * y2 - y3) / (6.0 * h);
}


double calc_diff2_m2_n4(double x2, double h) {
    double y0 = f(x2 - 2.0 * h);
    double y1 = f(x2 - h);
    double y2 = f(x2);
    double y3 = f(x2 + h);
    double y4 = f(x2 + 2.0 * h);
    return (-y0 + 16.0 * y1 - 30.0 * y2 + 16.0 * y3 - y4) / (12.0 * h * h);
}

class NumericalDifferentiator {
private:
    double x_min = -0.8;
    double x_max = 0.8;
    double step_h = 0.05;
    Gnuplot gp;

public:
    void set_h(double new_h) {
        step_h = new_h;
    }

    void plot_first_derivative() {
        std::vector<std::pair<double, double>> exact_pts;
        std::vector<std::pair<double, double>> num_pts;

        for (double x = x_min; x <= x_max; x += 0.005) {
            exact_pts.push_back({ x, df_exact(x) });
            num_pts.push_back({ x, calc_diff1_m1_n3(x, step_h) });
        }

        gp << "set title 'first derivative: m=1, n=3 in x_1 (h = " << step_h << ")'\n";
        gp << "set xlabel 'x'\n";
        gp << "set ylabel 'fd(x)'\n";
        gp << "set grid\n";
        gp << "plot '-' with lines title 'Precise fd(x)' lc rgb 'black' lw 2, "
            << "'-' with lines title 'Numeric (m=1, n=3)' lc rgb 'blue' dt 2 lw 2\n";

        gp.send1d(exact_pts);
        gp.send1d(num_pts);
        gp.flush();
        std::cout << "Press Enter to continue..." << std::endl;
        std::cin.get();
    }

    void plot_second_derivative() {
        std::vector<std::pair<double, double>> exact_pts;
        std::vector<std::pair<double, double>> num_pts;

        for (double x = x_min; x <= x_max; x += 0.005) {
            exact_pts.push_back({ x, d2f_exact(x) });
            num_pts.push_back({ x, calc_diff2_m2_n4(x, step_h) });
        }

        gp << "set title 'Second derivative: m=2, n=4 in x_2 (h = " << step_h << ")'\n";
        gp << "set xlabel 'x'\n";
        gp << "set ylabel 'fdd(x)'\n";
        gp << "set grid\n";
        gp << "plot '-' with lines title 'Precise fdd(x)' lc rgb 'black' lw 2, "
            << "'-' with lines title 'Numeric (m=2, n=4)' lc rgb 'red' dt 2 lw 2\n";

        gp.send1d(exact_pts);
        gp.send1d(num_pts);
        gp.flush();
        std::cout << "Press Enter to continue..." << std::endl;
        std::cin.get();
    }

    void print_error_table(double x_test) {
        std::cout << std::scientific << std::setprecision(6);
        std::cout << "\n=========================================================================\n";
        std::cout << "  Analysis of the error of numerical differentiation at point x = " << x_test << "\n";
        std::cout << "=========================================================================\n\n";

        std::cout << "1) First deirvative (m=1, n=3, theoretical order O(h^3)):\n";
        std::cout << "-------------------------------------------------------------------------\n";
        std::cout << "      Step h     |   Numeric   |    Precise     |  Deviation  | Ratio \n";
        std::cout << "-------------------------------------------------------------------------\n";

        double exact1 = df_exact(x_test);
        double prev_err1 = 0.0;
        double current_h = 0.1;

        for (int i = 0; i < 6; ++i) {
            double approx = calc_diff1_m1_n3(x_test, current_h);
            double err = std::abs(approx - exact1);
            std::cout << std::setw(14) << current_h << " | "
                << std::setw(13) << approx << " | "
                << std::setw(13) << exact1 << " | "
                << std::setw(13) << err << " | ";
            if (i > 0 && err > 0) {
                std::cout << std::fixed << std::setprecision(2) << std::setw(9) << (prev_err1 / err);
                std::cout << std::scientific << std::setprecision(6);
            }
            else {
                std::cout << std::setw(9) << "-";
            }
            std::cout << "\n";
            prev_err1 = err;
            current_h /= 2.0;
        }

        std::cout << "-------------------------------------------------------------------------\n\n";

        std::cout << "2) Second derivative (m=2, n=4, theoretical order O(h^4)):\n";
        std::cout << "-------------------------------------------------------------------------\n";
        std::cout << "      Step h     |   Numeric   |    Precise     |  Deviation  | Ratio \n";
        std::cout << "-------------------------------------------------------------------------\n";

        double exact2 = d2f_exact(x_test);
        double prev_err2 = 0.0;
        current_h = 0.1;

        for (int i = 0; i < 6; ++i) {
            double approx = calc_diff2_m2_n4(x_test, current_h);
            double err = std::abs(approx - exact2);
            std::cout << std::setw(14) << current_h << " | "
                << std::setw(13) << approx << " | "
                << std::setw(13) << exact2 << " | "
                << std::setw(13) << err << " | ";
            if (i > 0 && err > 0) {
                std::cout << std::fixed << std::setprecision(2) << std::setw(9) << (prev_err2 / err);
                std::cout << std::scientific << std::setprecision(6);
            }
            else {
                std::cout << std::setw(9) << "-";
            }
            std::cout << "\n";
            prev_err2 = err;
            current_h /= 2.0;
        }
        std::cout << "-------------------------------------------------------------------------\n\n";
    }
};

int main() {
    NumericalDifferentiator diff;

    diff.print_error_table(0.5);

    diff.set_h(0.05);
    diff.plot_first_derivative();
    diff.plot_second_derivative();

    return 0;
}