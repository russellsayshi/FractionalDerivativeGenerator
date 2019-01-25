#include <iostream>
#include <interpreter.h>
#include <lexer.h>
#include <iomanip>
#include <limits>
#include <cmath>

#define DELTA_X 0.0001
#define TAKE_AT_X DELTA_X * 100
#ifdef DEBUG
#define debug(x) std::cout << x << std::endl
#define ndebug(x) std::cout << x
#else
#define debug(x) do {} while(0)
#define ndebug(x) do {} while(0)
#endif

long double nth_derivative_of_monomial(long double power, long double deriv, long double position) {
	return std::tgamma(power+1)/std::tgamma(power-deriv+1)*pow(position, power-deriv);
}

int main(int argc, char** argv) {
	//Capture user data
	std::cout << "Enter a function of x: ";
	std::string func;
	getline(std::cin, func);
	long double nth_derivative = -1;
	while(nth_derivative <= 0 || nth_derivative > 1) {
		std::cout << "Enter the nth derivative you want to take (0 < n <= 1), or -1 to display Taylor series: ";
		std::cin >> nth_derivative;
		if(nth_derivative == -1) break;
	}
	std::cout << "Enter the x-coordinate at which to start: ";
	long double start_point;
	std::cin >> start_point;
	long num_taylor_terms = 3;
	std::cout << "Enter the number of points to output: ";
	int num_points_out;
	std::cin >> num_points_out;
	std::cout << "Enter the delta between each point: ";
	long double btw_point_delta;
	std::cin >> btw_point_delta;
	/*while(num_taylor_terms <= 0) {
		std::cout << "Enter the number of terms you want in the taylor series (recommendation ~80): ";
		std::cin >> num_taylor_terms;
	}*/

	//Create lexer and feed into interpreter
	lexer lex(func);
	interpreter inter;
	inter.fetch_tokens(lex);

	//Pick display settings
	std::cout << std::fixed;
	std::cout << std::setprecision(std::numeric_limits<long double>::digits10 + 1);

	//Assign variables
	std::unordered_map<std::string, long double> vars;
	vars["pi"] = 3.14159265359;
	vars["e"] = 2.71828182846;

	long double value_table[num_taylor_terms+1];
	long double derivatives[num_taylor_terms];

	for(int num_pts_computed = 0; num_pts_computed < num_points_out; num_pts_computed++) {
		long double point = start_point - (nth_derivative == -1 ? 0 : TAKE_AT_X) + num_pts_computed * btw_point_delta;
		for(int i = 0; i < num_taylor_terms+1; i++) {
			vars["x"] = point + DELTA_X * i;
			long double temp = inter.interpret(vars);
			if(std::isinf(temp) || std::isnan(temp)) {
				std::cerr << "Abort: at x = " << vars["x"] << " the value of f(x) was: " << temp << "." << std::endl;
				return 1;
			}
			value_table[i] = temp;
		}

		for(int i = 0; i < num_taylor_terms; i++) {
			derivatives[i] = value_table[0];

			ndebug(i << ": ");
			for(int o = 0; o < num_taylor_terms+1-i; o++) {
				ndebug(value_table[o] << "    ");
			}
			debug("");

			long double prev = value_table[0];
			for(int o = 0; o < num_taylor_terms - i; o++) {
				long double tmp = value_table[o+1];
				value_table[o] = (tmp - prev)/DELTA_X;
				prev = tmp;
			}
			//std::cout << "derivative[" << i << "] = " << derivatives[i] << std::endl;
		}

		if(nth_derivative == -1) {
			for(int i = 0; i < num_taylor_terms; i++) {
				long double coeff = 1/std::tgamma(i+1);
				std::cout << (derivatives[i] * coeff) << "*(x" << (point < 0 ? "+" : "-") << (point < 0 ? -point : point) << ")^" << i << (i == num_taylor_terms-1 ? "" : " + ");
			}
			std::cout << std::endl;
		} else {
			long double n_val = 0;
			for(int i = 0; i < num_taylor_terms; i++) {
				//long double t_coeff = 1 / std::tgamma(i+2);
				//long double nd_coeff = std::tgamma(i+2)/std::tgamma(i+2-nth_derivative);
				long double a_coeff = 1/std::tgamma(i+1-nth_derivative);
				long double x_term = pow(TAKE_AT_X, i-nth_derivative);
				debug("term " << i << ": " << derivatives[i] << "*" << a_coeff << "*x^(" << (i-nth_derivative) << ") = " << (a_coeff * x_term * derivatives[i]));
				n_val += a_coeff * x_term * derivatives[i];
				for(int o = -20; o < 20; o++) {
					std::cout << "(" << (o/4.0) << ", " << (pow(TAKE_AT_X + o/4.0, i-nth_derivative)*derivatives[i]*a_coeff) << "), ";
				}
				return 0;
			}
			std::cout << "(" << (point + TAKE_AT_X) << ", " << n_val << ")";
			if(num_pts_computed == num_points_out-1) std::cout << std::endl;
			else std::cout << ", ";
		}
	}
}
