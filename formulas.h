#pragma once
#include <cmath>
double a = 1.1, b = 2.3,
	   alpha = 4 / 5, beta = 0;

double f(double x)
{
	return 3.5 * cos(0.7 * x) * exp(-5.0 * x / 3.0)
		+ 2.4 * sin(5.5 * x) * exp(-3.0 * x / 4) + 5;
}

double leftRectangle(int nodes)
{
	double h = (b - a) / nodes;
	double integral_value = 0;
	for (int i = 0; i < (nodes - 1); ++i)
	{
		integral_value += f(a + h * i) * h;
	}
	return integral_value;
}

double middleRectangle(int nodes)
{
	double h = (b - a) / nodes;
	double integral_value = 0;
	for (int i = 0; i < (nodes - 1); ++i)
	{
		integral_value += f((2*a + h * (2*i + 1))/2) * h;
	}
	return integral_value;
}

double trapezoid(int nodes)
{
	double h = (b - a) / nodes;
	double integral_value = 0;
	for (int i = 0; i < (nodes - 1); ++i)
	{
		integral_value += (f(a + h*i) + f(a + h * (i+1))) * h / 2;
	}
	return integral_value;
}

double Simpson(int nodes)
{
	double h = (b - a) / nodes;
	double integral_value = 0;
	for (int i = 0; i < (nodes - 1); ++i)
	{
		integral_value += (f(a + h * i) 
			+ f(a + h * (i + 1)) 
			+ 4 * f((2 * a + h * (2 * i + 1)) / 2))
			* h / 6;
	}
	return integral_value;
}