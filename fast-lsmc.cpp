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
int num_trials, num_divisions, num_sims, call_flag;

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
Eigen::MatrixXd polynomial_regression(Eigen::MatrixXd independent, Eigen::MatrixXd dependent, const int order, const int num_obs, const string basis_type) {
    // Initialize matrices
	Eigen::MatrixXd X(num_obs, order);
	Eigen::MatrixXd y(num_obs, 1);

    // Fill in X based on the basis function.
	// This control flow is verbose but faster since it doesn't have to check basis_type multiple times.
	// See page 6 of https://jfin-swufe.springeropen.com/articles/10.1186/s40854-015-0019-0 for coefficient definitions.
	if (basis_type == "Power") {
		for (int j = 0; j < order; j++) {			// filling up columns
			for (int i = 0; i < num_obs; i++) {		// filling up rows
				X(i, j) = pow(independent(i, 0), j);
			}
    	}
	} else if (basis_type == "Laguerre") {
		for (int j = 0; j < order; j++) {
			for (int i = 0; i < num_obs; i++) {
				double eval = 0;
				for (int m = 0; m < j + 1; m++) {
					eval += pow(-1, m) * n_choose_k(j, m) / tgamma(m + 1) * pow(independent(i, 0), m);
				}
				X(i, j) = eval;
			}
    	}
	} else if (basis_type == "Legendre") {
		for (int j = 0; j < order; j++) {
			for (int i = 0; i < num_obs; i++) {
				double eval = 0;
				for (int m = 0; m < (j / 2) + 1; m++) {
					eval += pow(-1, m) * n_choose_k(j, m) * n_choose_k(static_cast<int>(2 * (j - m)), j) * pow(independent(i, 0), j - 2 * m);
				}
				X(i, j) = eval * pow(1 / 2, j);
			}
    	}
	} else if (basis_type == "Hermitian") {
		for (int j = 0; j < order; j++) {
			for (int i = 0; i < num_obs; i++) {
				double eval = 0;
				for (int m = 0; m < (j / 2) + 1; m++) {
					eval += pow(-1, m) * pow(2 * independent(i, 0), j - 2 * m) / (tgamma(m + 1) * tgamma(static_cast<int>(j - 2 * m) + 1));
				}
				X(i, j) = eval * tgamma(j + 1);
			}
    	}
	}
    
	
	// Fill in y
	for (int i = 0; i < num_obs; i++) {
        y(i, 0) = dependent(i, 0);
    }
	
    // Solving
	return (X.transpose() * X).inverse() * X.transpose() * y;
}

double adjusted_determination_coef(const double observations[], const double predictions[], const int num_variables) {
	int num_observations = sizeof(observations) / sizeof(observations[0]);

	// Calculating observation mean for SST
	double observation_mean = 0;
	for (int i = 0; i < num_observations; i++) {
		observation_mean += observations[i];
	}
	observation_mean /= num_observations;

	// Calculating Total Sum of Squares (SST)
	double SST = 0;
	for (int i = 0; i < num_observations; i++) {
		SST += pow(observations[i] - observation_mean, 2);
	}

	// Calculating Residual Sum of Squares (RSS)
	double SSE = 0;
	for (int i = 0; i < num_observations; i++) {
		SSE += pow(observations[i] - predictions[i], 2);
	}

	// R^2_{adj} formula
	return 1 - (1 - SSE / SST) * ((num_observations - 1) / (num_observations - 1 - num_variables));
}

// n_C_k function: https://stackoverflow.com/questions/9330915/number-of-combinations-n-choose-r-in-c
int n_choose_k(const int n, int k) {
    if (k > n) {
		return 0;
	}
    if (k * 2 > n) {
		k = n - k;
	}
    if (k == 0) {
		return 1;
	}

    int result = n;
    for(int i = 2; i <= k; i++) {
        result *= (n - i + 1);
        result /= i;
    }
    return result;
}

// See page 6 of https://jfin-swufe.springeropen.com/articles/10.1186/s40854-015-0019-0 for coefficient definitions.
// Legendre coefficient function. Defined on [-1, 1].
// double legendre_coeffs(const int index, const int max_index) {
// 	return pow(-2, max_index) * pow(-1, index) * n_choose_k(max_index, index) * n_choose_k(static_cast<int>(2 * (index - max_index)), max_index);
// }

// // Laguerre polynomial coeffient function. Defined on [0, inf].
// double laguerre_coeffs(const int index, const int max_index) {
// 	return n_choose_k(max_index, index) * pow(-1, index) / tgamma(index + 1); // uses gamma function + 1, equiv to factorial
// }

// // Hermitian polynomial coeffcient function. Defined on [-inf, inf].
// double hermitian_coeffs(const int index, const int max_index) {
// 	return tgamma(max_index + 1) * pow(-1, index) / (tgamma(index + 1) * tgamma(static_cast<int>(max_index - 2 * index + 1)));
// }

int main(int argc, char* argv[]) {
	// reading arguments
    sscanf(argv[1], "%lf", &T); 			// expiry time
	sscanf(argv[2], "%lf", &r); 			// risk-free rate
	sscanf(argv[3], "%lf", &sigma); 		// volatility
	sscanf(argv[4], "%lf", &S_0);			// initial price
	sscanf(argv[5], "%lf", &K);				// strike price
	sscanf(argv[6], "%d", &num_divisions);	// # grid steps
	sscanf(argv[7], "%d", &num_sims);		// number simulations
	sscanf(argv[8], "%d", &call_flag);		// call (1) or put (0)

    // defining time, rate, and SD step lengths
    double dt = T / ((double) num_divisions);
	double dR = (r - 0.5 * pow(sigma, 2)) * dt;
	double dSD = sigma * sqrt(dt);
	double R = exp(r * T / ((double) num_divisions));

    // printing arguments (sanity check)
    cout << "--------------------------------" << endl;
	cout << "American Option Pricing using (Modified) Longstaff and Schwartz's Least Squares Monte Carlo Simulation" << endl;
	cout << "Expiration Time (Years) = " << T << endl;
	cout << "Risk Free Interest Rate = " << r << endl;
	cout << "Volatility (%age of stock value) = " << sigma * 100 << endl;
	cout << "Initial Stock Price = " << S_0 << endl;
	cout << "Strike Price = " << K << endl;
	cout << "Number of Simulations = " << num_sims << endl;
	cout << "Number of Divisions = " << num_divisions << endl;
	cout << "R = " << R << endl;
	cout << "Is Call? = " << call_flag << endl;
	cout << "--------------------------------" << endl;
    
    // defining American put price
    double option_price = 0.0;

    // Iterating through blocks, this is to limit memory usage since solving using num_sims would use big matrices.
	// The outcome is the same doing it this way.
    for (int k = 0; k <= (num_sims / 200); k++) {
		
		// setting trial number
		if (k != num_sims / 200) {
            num_trials = 200;
        } else {
            num_trials = num_sims % 200;
        }			
		
		if (num_trials > 0) {
			double asset_price[num_trials][num_divisions];

			// Setting initial price
			for (int i = 0; i < num_trials; i++) {
                asset_price[i][0] = S_0;
            }
			
			// Generating random path
			for (int i = 0; i < num_trials; i++) {
                for (int j = 1; j < num_divisions; j++) {
                    asset_price[i][j] = asset_price[i][j - 1] * exp(dR + dSD * get_gaussian());
                }
            }
            
			double value[num_trials];

			// Setting last column to payoff instead of price
			for (int i = 0; i < num_trials; i++) {
				if (call_flag == 1) {
					value[i] = max(0.0, asset_price[i][num_divisions - 1] - K);
				} else {
					value[i] = max(0.0, K - asset_price[i][num_divisions - 1]);
				}
            }
			
			for (int i = (num_divisions - 1); i > 0; i--) {
				Eigen::MatrixXd independent_vars(num_trials, 1);
				Eigen::MatrixXd dependent_vars(num_trials, 1);

				int num_paths = 0;

				// Choosing which values to use for regression/interpolation
				// Adapting this to non-zero cashflow. If current rebate + discounted cashflow up to exercise < discounted min rebate
				// at next exercise date, then do not evaluate. This is changed from if K - S_i > 0 then keep.
				// This is basically just another method of sub-optimal point elimination that provides a greater lower bound.

				// Taken from last sentence of first paragraph of page 12 of:
				// https://papers.ssrn.com/sol3/papers.cfm?abstract_id=1331904
				for (int j = 0; j < num_trials; j++) {
					if (call_flag == 1) {
						if ((asset_price[j][i] - K) + (asset_price[j][i + 1] - K) / R * exp(r * i / num_divisions * T) < 
							(asset_price[j][i + 1] - K) / R * exp(r * (i + 1) / num_divisions * T)) {
							num_paths++;
							independent_vars(num_paths, 0) = asset_price[j][i];
							dependent_vars(num_paths, 0) = value[j] / R;
						}
					} else {
						if ((K - asset_price[j][i]) + (K - asset_price[j][i + 1]) / R * exp(r * i / num_divisions * T) > 
							(K - asset_price[j][i + 1]) / R * exp(r * (i + 1) / num_divisions * T)) {
							num_paths++;
							independent_vars(num_paths, 0) = asset_price[j][i];
							dependent_vars(num_paths, 0) = value[j] / R;
						}
					}
				}
				
				if (num_paths > 0) {
					// regressing the dependent_variables on the independent variables using a 4th order polynomial
					// Power basis; change this later. This provides a lower bound (for call) but not an upper bound (for call)

					// The bases that will be checked are Power, Legendre, Laguerre, and Hermite A.
					int poly_degree = min(num_paths, 10);
					Eigen::MatrixXd a(poly_degree, 1);
					a = polynomial_regression(independent_vars, dependent_vars, poly_degree, num_paths, "Power");
					double greatest_r_sq_adj = 0;
					double optimal_poly_eval;
					
					for (int j = 0; j < num_trials; j++) {
						double power_poly_eval = 0;
						// double legendre_poly_eval = 0;
						// double laguerre_poly_eval = 0;
						// double hermite_poly_eval = 0;

						// Calculating polynomial evaluations of each.
						// Each path will choose the optimal basis by maximizing R^2_adj.

						// THE BASIS NEEDS TO BE TAKEN CARE OF IN THE POLYNOMIAL REGRESSION FUNCTION
						for (int l = 0; l < poly_degree; l++) {
							power_poly_eval += a(l, 0) * pow(asset_price[j][i], l);
							// laguerre_poly_eval += laguerre_coeffs(l, poly_degree - 1) * pow(asset_price[j][i], l);
							// if (l <= poly_degree / 2) {
							// 	legendre_poly_eval += legendre_coeffs(l, poly_degree - 1) * pow(asset_price[j][i], poly_degree - 1 - 2 * l);
							// 	hermite_poly_eval += hermitian_coeffs(l, poly_degree - 1) * pow(2 * asset_price[j][i], poly_degree - 1 - 2 * l);
							// }
						}
						
						adjusted_determination_coef;
						optimal_poly_eval = 4;

						if (call_flag == 1) {
							// (S - K) > poly_eval is the Andersen trigger method for which convergence is sped up.
							// This maximizes the average cutoff over the simulated paths.
							if (((asset_price[j][i] - K) > optimal_poly_eval) && ((asset_price[j][i] - K) > 0.0)) {
								value[j] = asset_price[j][i] - K;
							} else {
								value[j] /= R;
							}
						} else {
							if (((K - asset_price[j][i]) > optimal_poly_eval) && ((K - asset_price[j][i]) > 0.0)) {
								value[j] = K - asset_price[j][i];
							} else {
								value[j] /= R;
							}
						}
					}
				}
			}
			
			double local_option_price = 0.0;

			for (int j = 0; j < num_trials; j++) {
                local_option_price += value[j];
            }

			local_option_price /= ((double) num_trials) * R;
			option_price += local_option_price;
		}
	}
	
	if (num_sims % 200 == 0) 
		cout << "American Option Price = " << option_price / ((double) num_sims / 200)  << endl;
	else 
		cout << "Option Price = " << option_price / ((double) num_sims / 200 + 1) << endl;
}
// Hypothetically, this should be more accurate and faster than the standard LS. I'm a complete imbecile so maybe not.