#pragma once
#include <vector>
std::vector<double> gauss(std::vector<std::vector<double>> matrixx, std::vector<double> col, int n) // метод Гаусса
{
    std::vector<double> result(n, 0);
    int i, j, l;
    double boo;
    for (l = 0; l < n; l++) // прямой ход
    {
        if (matrixx[l][l] == 0)
        {
            for (int t = l + 1; t < n; t++)
            {
                if (matrixx[t][l] != 0) std::swap(matrixx[t], matrixx[l]);
            }
        }
        for (i = l + 1; i < n; i++)
        {
            boo = matrixx[i][l];
            for (j = l; j < n; j++)
            {
                matrixx[i][j] -= (matrixx[l][j] / matrixx[l][l]) * boo;
            }
            col[i] -= (col[l] / matrixx[l][l]) * boo;
        }
    }
    result[n - 1] = col[n - 1] / matrixx[n - 1][n - 1];

    for (i = n - 2; i >= 0; i--) // обратный
    {
        result[i] = col[i];
        for (j = n - 1; j > i; j--)
        {
            result[i] -= matrixx[i][j] * result[j];
        }
        result[i] /= matrixx[i][i];
    }
    return result;
}