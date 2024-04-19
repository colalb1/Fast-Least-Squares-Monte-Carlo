// Importing packages
#include <iostream>
#include <iomanip>

#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <random>
#include <chrono>
#include <eigen3/Eigen/Dense>

// Declaring namespaces
using namespace Eigen;
using namespace std;

// Setting random seed
unsigned seed = (unsigned) std::chrono::system_clock::now().time_since_epoch().count();
default_random_engine generator;

// Initializing variables
double r, K, S_0, T, sigma;
int num_trials, num_divisions, num_sims;

// Max and min functions
double max(double a, double b) {
	return (b < a) ? a:b;
}

int min(int a, int b) {
	return (b < a) ? b:a;
}

// Generating standard uniform
double get_uniform() {
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    return distribution(generator);
}

// Box-Muller transform
double get_gaussian() {
	return sqrt(-2.0 * log(get_uniform())) * cos(2 * M_PI * get_uniform());
}

// Polynomial regression using (X^T * X)^{-1} * X^T * y
Eigen::MatrixXd polynomial_regression(Eigen::MatrixXd independent, Eigen::MatrixXd dependent, int order, int num_obs) {
    // initialize matrices
	Eigen::MatrixXd X(num_obs, order);
	Eigen::MatrixXd y(num_obs, 1);

    // fill in X and y
    for (int j = 0; j < order; j++) {
        for (int i = 0; i < num_obs; i++) {
            X(i, j) = pow(independent(i, 0), j - 1);
        }
    }
	
	for (int i = 0; i < num_obs; i++) {
        y(i, 0) = dependent(i, 0);
    }
	
    // solving
	return (X.transpose() * X).inverse() * X.transpose() * y;
}

int main(int argc, char* argv[]) {
	// reading arguments
    sscanf(argv[1], "%lf", &T); 			// expiry time
	sscanf(argv[2], "%lf", &r); 			// risk-free rate
	sscanf(argv[3], "%lf", &sigma); 		// volatility
	sscanf(argv[4], "%lf", &S_0);			// initial price
	sscanf(argv[5], "%lf", &K);				// strike price
	sscanf(argv[6], "%d", &num_divisions);	// # grid steps
	sscanf(argv[7], "%d", &num_sims);		// number simulations

    // defining time, rate, and SD step lengths
    double dt = T / ((double) num_divisions);
	double dR = (r - 0.5 * pow(sigma, 2)) * dt;
	double dSD = sigma * sqrt(dt);
	double R = exp(r * T / ((double) num_divisions));

    // printing arguments (sanity check)
    cout << "--------------------------------" << endl;
	cout << "American Put Option Pricing using Longstaff and Schwartz's Least Squares Monte Carlo Simulation" << endl;
	cout << "Expiration Time (Years) = " << T << endl;
	cout << "Risk Free Interest Rate = " << r << endl;
	cout << "Volatility (%age of stock value) = " << sigma * 100 << endl;
	cout << "Initial Stock Price = " << S_0 << endl;
	cout << "Strike Price = " << K << endl;
	cout << "Number of Simulations = " << num_sims << endl;
	cout << "Number of Divisions = " << num_divisions << endl;
	cout << "R = " << R << endl;
	cout << "--------------------------------" << endl;
    
    // defining American put price
    double put_price = 0.0;

    // iterating through blocks, this is to limit memory usage since solving using num_sims would use big matrices
	// the outcome is the same doing it this way
    for (int k = 0; k <= (num_sims / 200); k++) {
		
		// setting trial number
		if (k != num_sims / 200) {
            num_trials = 200;
        } else {
            num_trials = num_sims % 200;
        }			
		
		if (num_trials > 0) {
			double asset_price[num_trials][num_divisions];

			for (int i = 0; i < num_trials; i++) {
                asset_price[i][0] = S_0;
            }
			
			for (int i = 0; i < num_trials; i++) {
                for (int j = 1; j < num_divisions; j++) {
                    asset_price[i][j] = asset_price[i][j - 1] * exp(dR + dSD * get_gaussian());
                }
            }
            
			double value[num_trials];

			// initialize the value based on the price at final stage
			for (int i = 0; i < num_trials; i++) {
                value[i] = max(0.0, K - asset_price[i][num_divisions - 1]);
            }
			
			for (int i = (num_divisions - 1); i > 0; i--) {
				Eigen::MatrixXd independent_vars(num_trials, 1);
				Eigen::MatrixXd dependent_vars(num_trials, 1);

				int num_vars = 0;

				for (int j = 0; j < num_trials; j++) {
					if (max(0.0, K - asset_price[j][i]) > 0) {
						num_vars++;
						independent_vars(num_vars, 0) = asset_price[j][i];
						dependent_vars(num_vars, 0) = value[j] / R;
					}
				}
				
				if (num_vars > 0) {
					// regressing the dependent_variables on the independent variables using a 4th order polynomial
					Eigen::MatrixXd a(min(num_vars, 5), 1);
					a = polynomial_regression(independent_vars, dependent_vars, min(5, num_vars), num_vars);
                    int poly_degree = min(num_vars, 5);

					for (int j = 0; j < num_trials; j++) {
						double poly_eval = 0;

						for (int l = 0; l < poly_degree; l++) {
							poly_eval += a(l, 0) * pow(asset_price[j][i], l);
						}
						
						if (((K - asset_price[j][i]) > poly_eval) && ((K - asset_price[j][i]) > 0.0)) {
							value[j] = K - asset_price[j][i];
						} else {
							value[j] /= R;
						}
					}
				}
			}
			
			double local_put_price = 0.0;

			for (int j = 0; j < num_trials; j++) {
                local_put_price += value[j];
            }

			local_put_price = (local_put_price / ((double) num_trials)) / R;
			put_price += local_put_price;
		}
	}
	
	if (num_sims % 200 == 0) 
		cout << "American Put Price = " << put_price / ((double) num_sims / 200)  << endl;
	else 
		cout << "Put Price = " << put_price / ((double) num_sims / 200 + 1) << endl;
}