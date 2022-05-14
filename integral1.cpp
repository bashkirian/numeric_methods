#include "formulas.h"
#include "newton_gauss.h"
#include "iomanip"
#include <iostream>
#include <string>
#include <limits>
using namespace std;

const double EPS = pow(10, -6);
int nodes_error = 2;
double h_initial = 1.2;
vector<vector<double>> make_matrix(int r, int m, double h)
{
    vector<vector<double>> result; 
    for (int i = 0; i < (r + 1); ++i)
    {
        vector<double> row;
        for (int j = 0; j < r; ++j)
        {
            row.push_back(-pow((h_initial / pow(2, i)), m + j));
        }
        row.push_back(1.0);
        result.push_back(row);
    }
    return result;
}

void SKFError(SKF& skf, double h, string& method)
{
    double error = numeric_limits<double>::max();
    vector<double> integrals;
    vector<double> coefficients;
    vector<vector<double>> matrix;
    skf.SetStep(h);
    for (int i = 0; i < 3; ++i)
    {
        if (method == "Newton")
            integrals.push_back(skf.Newton_Cotes());
        else
            integrals.push_back(skf.Gauss());
        skf.SetStep(h / pow(2, (i + 1)));
    }
    int m = -log(abs((integrals[2] - integrals[1]) / (integrals[1] - integrals[0]))) / log(2);
    int r = 1;
	double h_opt = h * pow((EPS * (1 - pow(2, -m))) / abs(integrals[1] - integrals[0]), 1 / m);
	skf.SetStep(h_opt);
    integrals.pop_back();
    cout << setw(15) << "Шаг" << setw(25) << setprecision(10) << "Вычисленная кв.сумма" << setw(20) 
        << "Погрешность" 
        << setw(10) << "Порядок" << endl;
    double integral_value;
    while (true)
    {
        matrix = make_matrix(r, m, h);
        vector<double> coefficients = gauss(matrix, integrals, r + 1);
        error = abs(coefficients.back() - integrals[r - 1]);
        cout << setw(15) << h << setw(25) << setprecision(10) << coefficients.back() << setw(20) << error;
        if (error < EPS) break;
        if (method == "Newton")
            integrals.push_back(skf.Newton_Cotes());
        else
            integrals.push_back(skf.Gauss());
        nodes_error /= 2;
        h /= 2;
        cout << setw(10) << r << endl;
        r++;
    }
    cout << "Оптимальный шаг = " << h << endl;
}

int main()
{
    int nodes;
    setlocale(LC_ALL, "Russian");
    SKF skf(1);
    cout << setw(15) << "Кол-во сегментов"
        << setw(15) << "Левого прям-ка"
        << setw(15) << "Среднего прям-ка"
        << setw(15) << "Трапеции"
        << setw(15) << "Симпсона"
        << setw(15) << "Ньютона-Котеса"
        << setw(15) << "Гаусса" << endl;
    for (int segments = 5; segments <= 105; segments += 10)
    {
        skf.SetSegments(segments);
        cout << setw(15) << segments;
        cout << setw(15) << leftRectangle(segments);
        cout << setw(15) << middleRectangle(segments);
        cout << setw(15) << trapezoid(segments);
        cout << setw(15) << Simpson(segments);
        cout << setw(15) << skf.Newton_Cotes();
        cout << setw(15) << skf.Gauss() << endl;
    }
    double h = h_initial;
    cout << endl << endl << "Оптимальный шаг для Ньютона-Котеса с начальным шагом " << h << endl;
    string method = "Newton";
    SKFError(skf, h, method);
    cout << endl << endl << "Оптимальный шаг для Гаусса с начальным шагом " << h << endl;
    method = "Gauss";
    SKFError(skf, h, method);
}