#pragma once
#include "SLAE.h"
#include <cmath>
#include <vector>

class SKF
{
	double a = 1.1, b = 2.3,
		alpha = 4 / 5;
	std::vector<double> nodes;
	double f(double x)
	{
		return 3.5 * cos(0.7 * x) * exp(-5.0 * x / 3.0)
			+ 2.4 * sin(5.5 * x) * exp(-3.0 * x / 4) + 5;
	}

	double moment(size_t i, size_t order)
	{
 		double result = pow((nodes[i] - a), order + 1 - alpha)
			- pow((nodes[i - 1] - a), order + 1 - alpha);
		result /= (order + 1 - alpha);
		switch (order)
		{
		case 0:
			return result;
		case 1:
			return result + a * moment(i, 0);
		case 2:
			return result + 2 * a * moment(i, 1)
				- a * a * moment(i, 0);
		case 3:
			return result + 3 * a * moment(i, 2) - 3 * pow(a, 2) * moment(i, 1) + pow(a,3) * moment(i, 0);
		case 4:
			return result + 4 * a * moment(i, 3)
				- 6 * pow(a, 2) * moment(i, 2) +
				4 * pow(a, 3) * moment(i, 1)
				- pow(a, 4) * moment(i, 0);
		case 5:
			return result + 5 * a * moment(i, 4)
				- 10 * pow(a, 2) * moment(i, 3)
				+ 10 * pow(a, 3) * moment(i, 2)
				- 5 * pow(a, 4) * moment(i, 1)
				+ pow(a, 5) * moment(i, 0);
		}
	}

	std::vector<double> cardano(const std::vector<double>& coeffs)
	{
		std::vector<double> result(3);
		double a_ = 1;
		double b_ = coeffs[2];
		double c = coeffs[1];
		double d = coeffs[0];
		double q = 2 * pow(b_, 3) / (27 * pow(a_, 3)) - b_ * c / (3 * a_ * a_) + d / a_;
		q /= 2;
		double p = (3 * a_ * c - b_ * b_) / (3 * a_ * a_) / 3;
		double r = sqrt(abs(p)) * q / abs(q);
		double phi = acos(q / pow(r, 3));
		result[0] = -2 * r * cos(phi / 3) - b_ / 3;
		result[1] = 2 * r * cos(3.14 / 3 - phi / 3) - b_ / 3;
		result[2] = 2 * r * cos(3.14 / 3 + phi / 3) - b_ / 3;
		return result;
	}

public:
	SKF(int number_of_segments)
	{
		double h = (b - a) / number_of_segments;
		for (int i = 0; i < (number_of_segments + 1); ++i)
		{
			nodes.push_back(a + h * i);
		}
	}

	void SetSegments(int number_of_segments)
	{
		double h = (b - a) / number_of_segments;
		nodes.clear();
		for (int i = 0; i < (number_of_segments + 1); ++i)
		{
			nodes.push_back(a + h * i);
		}
	}

	void SetStep(double h)
	{
		double number_of_segments = (b - a) / h;
		nodes.clear();
		for (int i = 0; i < (number_of_segments + 1); ++i)
		{
			nodes.push_back(a + h * i);
		}
	}

	double Newton_Cotes()
	{
		std::vector<double> moments(3);
		double integral = 0;
		std::vector<std::vector<double>> A_matrix(3, std::vector<double>(3, 0));
		std::vector<double> column(3, 0);
		std::vector<double> x_nodes(3);
		std::vector<double> A_coeffs(3);
		for (int i = 1; i < nodes.size(); ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				moments[j] = moment(i, j);
			}
			x_nodes[0] = nodes[i - 1];
			x_nodes[1] = (nodes[i - 1] + nodes[i]) / 2;
			x_nodes[2] = nodes[i];
			for (int row = 0; row < 3; ++row)
			{
				for (int col = 0; col < 3; ++col)
				{
					A_matrix[row][col] = pow(x_nodes[col], row);
				}
			}
			column[0] = moments[0];
			column[1] = moments[1];
			column[2] = moments[2];
			A_coeffs = gauss(A_matrix, column, 3);
			integral += A_coeffs[0] * f(x_nodes[0])
				+ A_coeffs[1] * f(x_nodes[1])
				+ A_coeffs[2] * f(x_nodes[2]);
		}
		return integral;
	}

	double Gauss()
	{
		std::vector<double> moments(6);
		double integral = 0;
		std::vector<std::vector<double>> A_matrix(3, std::vector<double>(3, 0));
		std::vector<double> column(3, 0);
		std::vector<double> result_a;
		std::vector<double> x_nodes(3);
		std::vector<double> A_coeffs(3);
		for (int i = 1; i < nodes.size(); ++i)
		{
			for (int j = 0; j <= 5; ++j)
			{
				moments[j] = moment(i, j);
			}
			for (int row = 0; row < 3; ++row)
			{
				for (int col = 0; col < 3; ++col)
				{
					A_matrix[row][col] = moments[row + col];
				}
			}
			column[0] = -moments[3];
			column[1] = -moments[4];
			column[2] = -moments[5];
			result_a = gauss(A_matrix, column, 3);
			x_nodes = cardano(result_a);
			for (int row = 0; row < 3; ++row)
			{
				for (int col = 0; col < 3; ++col)
				{
					A_matrix[row][col] = pow(x_nodes[col], row);
				}
			}
			column[0] = moments[0];
			column[1] = moments[1];
			column[2] = moments[2];
			A_coeffs = gauss(A_matrix, column, 3);
			integral += A_coeffs[0] * f(x_nodes[0])
				      + A_coeffs[1] * f(x_nodes[1])
				      + A_coeffs[2] * f(x_nodes[2]);
		}
		return integral;
	}
};
