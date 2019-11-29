#include <iostream>
#include <cmath>

using namespace std;
const double lambda = 100.0;
const double a = 0;
const double b = 1;

double f(double x, double u)
{
	return -lambda*u;
}
double U_exact(double x)
{
	return exp(-lambda*x);
}	// yold		//x	//z	
double F(double y0, double x1, double y1, double h )
{
	return y1 - y0 - h*f(x1, y1);
}
double F1(double h)
{
	// Generalize for x
	return 1 + lambda*h;	
}
double euler_explicit(double y0, double h);
double euler_implicit(double y0, double h);
double newton_method(double y0, double x, double h);
int main()
{
	double h = 0.1,  y0 =1;		
	cout << "U(" << b << ") =  " << U_exact(b) << endl;
	cout << euler_explicit(y0, h);
	cout << endl << endl;
	cout << "U(" << b << ") =  " << U_exact(b) << endl;
	cout << euler_implicit(y0, h);
	return 0;
}

double euler_explicit(double y0,  double h)
{
	double y_old = y0, y_new, x = a, eps=0.0001;
	do{
		std::cout<< "Approximation: at " << x << " " 
			<< y_old << std::endl;
		
		x = x + h;
	 	y_new = y_old + h* f(x, y_old);
		y_old = y_new;
		
	} while(fabs(b - x) > eps);
	return y_new;
}


double euler_implicit(double y0, double h)
{
	double y_new = y0, y_old = y0, x = a,  eps = 0.001;
	do{
		std::cout<< "Approximation: at " << x << " " 
			<< y_old << std::endl;
		x = x + h;
		y_old = y_new;
		y_new = newton_method(y_old, x, h);
	}while(fabs(b-x) > eps);
	return y_new;
}
double newton_method(double y0, double x, double h)
{
	int k = 0;
	double znew = y0,  zold ;
	do{
		k++;
		zold = znew;
		znew = zold - F(y0, x, zold, h)/F1(h);

	}while(fabs(znew - zold) > 0.0001 & k < 1500);
	if(k == 1500) cout << " k = 1500" << endl;	
	return znew;
}

