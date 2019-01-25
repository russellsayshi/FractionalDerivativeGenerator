#include <iostream>
#include <interpreter.h>
#include <lexer.h>
#include <cmath>

#define DELTA_X 0.0001
#define TAKE_AT_X = DELTA_X * 1000

int main(int argc, char** argv) {
	//Capture user data
	std::cout << "Enter a function of x: ";
	std::string func;
	getline(std::cin, func);
	long double nth_derivative = -1;
	while(nth_derivative <= 0 || nth_derivative > 1) {
		std::cout << "Enter the nth derivative you want to take (0 < n <= 1): ";
		std::cin >> nth_derivative;
	}
	std::cout << "Enter the x-coordinate at which to center the taylor series: ";
	long double point;
	std::cin >> point;
	long num_taylor_terms = 3;
	/*while(num_taylor_terms <= 0) {
		std::cout << "Enter the number of terms you want in the taylor series (recommendation ~80): ";
		std::cin >> num_taylor_terms;
	}*/

	//Create lexer and feed into interpreter
	lexer lex(func);
	interpreter inter;
	inter.fetch_tokens(lex);

	//Assign variables
	std::unordered_map<std::string, long double> vars;
	vars["pi"] = 3.14159265359;
	vars["e"] = 2.71828182846;

	long double value_table[num_taylor_terms+1];
	bool has_put_initial = false;
	long double initial;
	for(int i = 0; i < num_taylor_terms+1; i++) {
		vars["x"] = point + DELTA_X * i;
		long double temp = inter.interpret(vars);
		if(!has_put_initial) {
			initial = temp;
			has_put_initial = true;
		}
		value_table[i] = temp;
		if(std::isinf(temp) || std::isnan(temp)) {
			std::cerr << "Abort: at x = " << vars["x"] << " the value of f(x) was: " << temp << "." << std::endl;
			return 1;
		}
	}

	long double values_copy[num_taylor_terms+1];
	long double derivatives[num_taylor_terms];
	for(int i = 0; i < num_taylor_terms; i++) {
		derivatives[i] = value_table[0];
			std::cout << i << ": ";
			for(int o = 0; o < num_taylor_terms+1-i; o++) {
				std::cout << value_table[o] << "    ";
			}
			std::cout << std::endl;
		long double prev = value_table[0];
		for(int o = 0; o < num_taylor_terms - i; o++) {
			long double tmp = value_table[o+1];
			value_table[o] = (tmp - prev)/DELTA_X;
			prev = tmp;
		}
		//std::cout << "derivative[" << i << "] = " << derivatives[i] << std::endl;
	}

	double take_at_x = 0.1;
	double n_val = initial/std::tgamma(2-nth_derivative)*pow(take_at_x-point, -nth_derivative);
	for(int i = 0; i < num_taylor_terms; i++) {
		//long double t_coeff = 1 / std::tgamma(i+2);
		//long double nd_coeff = std::tgamma(i+2)/std::tgamma(i+2-nth_derivative);
		long double a_coeff = 1/std::tgamma(i+2-nth_derivative);
		long double x_term = pow(take_at_x - point, i+1-nth_derivative);
		n_val += a_coeff * x_term;
	}
	std::cout << n_val << std::endl;
}
