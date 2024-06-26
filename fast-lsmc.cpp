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
string basis;

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

// Polynomial regression using (X^T * X)^{-1} * X^T * y
tuple<Eigen::MatrixXd, Eigen::MatrixXd> polynomial_regression(Eigen::MatrixXd independent, Eigen::MatrixXd dependent, int order, int num_obs, string basis_type) {
    // Initialize matrices
	Eigen::MatrixXd X(num_obs, order);
	Eigen::MatrixXd y(num_obs, 1);

    // Fill in X based on the basis function.
	// This control flow is verbose but faster since it doesn't have to check basis_type multiple times.
	// See page 6 of https://jfin-swufe.springeropen.com/articles/10.1186/s40854-015-0019-0 for coefficient definitions.
	if (basis_type == "Power") {
		for (int j = 0; j < order; j++) {			// filling up columns
			for (int i = 0; i < num_obs; i++) {		// filling up rows
				X(i, j) = pow(independent(i, 0), j - 1);
			}
    	}
	} else if (basis_type == "Laguerre") {
		for (int j = 0; j < order; j++) {
			for (int i = 0; i < num_obs; i++) {
				double poly_eval = 0;

				for (int m = 0; m < j + 1; m++) {
					poly_eval += n_choose_k(j, m) / tgamma(m + 1) * pow(-independent(i, 0), m);
				}
				X(i, j) = poly_eval;
			}
    	}
	} else if (basis_type == "Hermitian") {
		for (int j = 0; j < order; j++) {
			for (int i = 0; i < num_obs; i++) {
				double poly_eval = 0;

				for (int m = 0; m < (j / 2) + 1; m++) {
					poly_eval += pow(-1, m) * pow(2 * independent(i, 0), j - 2 * m) / (tgamma(m + 1) * tgamma(static_cast<int>(j - 2 * m) + 1));
				}
				X(i, j) = tgamma(j + 1) * poly_eval;
			}
    	}
	}
    
	// Fill in y
	for (int i = 0; i < num_obs; i++) {
        y(i, 0) = dependent(i, 0);
    }
	
    // Solving and returning (X, beta)
	return make_tuple(X, (X.transpose() * X).inverse() * X.transpose() * y);
}

// This function is dead; was supposed to be used for selecting the basis programmatically but that was slow 
double adjusted_determination_coef(Eigen::MatrixXd observations, Eigen::MatrixXd predictions, int num_variables) {
	int num_observations = observations.rows();
	
	// Calculating observation mean for SST
	double observation_mean = 0;
	for (int i = 0; i < num_observations; i++) {
		observation_mean += observations(i, 0);
	}
	observation_mean /= num_observations;
	
	// Calculating Total Sum of Squares (SST)
	double SST = 0;
	for (int i = 0; i < num_observations; i++) {
		SST += pow(observations(i, 0) - observation_mean, 2);
	}
	
	// Calculating Residual Sum of Squares (RSS)
	double SSE = 0;
	for (int i = 0; i < num_observations; i++) {
		SSE += pow(observations(i, 0) - predictions(i, 0), 2);
	}
	
	// Interpolation case
	if ((num_observations - 1 - num_variables) == 0) {
		return 1.0;
	}
	// R^2_{adj} formula
	return 1 - (1 - SSE / SST) * ((num_observations - 1) / (num_observations - 1 - num_variables));
}

// Hermitian or Laguerre when T > 1 and K < S_0 overestimate the option price
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
	basis = argv[9];						// Power, Laguerre, Hermitian

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
	cout << "Volatility = " << sigma << endl;
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

	// Stratification and double-regression enhancement from https://www.sciencedirect.com/science/article/pii/S0165188913000493
	// cannot be done since each region is only one timestep long since this is an American option as opposed to 
	// a Bermudan option. This means I must implement other optimizations. Batched iteration also prevents this since I need
	// access to all of the paths.
    for (int k = 0; k <= (num_sims / 200); k++) {
		
		// setting trial number
		if (k != num_sims / 200) {
            num_trials = 200;
        } else {
            num_trials = num_sims % 200;
        }			
		
		if (num_trials > 0) {
			double brownian_bridge_values[num_trials], next_brownian_bridge_values[num_trials], value[num_trials];

            // Calculating last point
            for (int i = 0; i < num_trials; i++) {
                // Initializing end price
                brownian_bridge_values[i] = (r - pow(sigma, 2) / 2) * T + sigma * sqrt(T) * get_gaussian();

                // Calculating payoff
                if (call_flag == 1) {
                    value[i] = max(0.0, S_0 * exp(brownian_bridge_values[i]) - K);
                } else {
                    value[i] = max(0.0, K - S_0 * exp(brownian_bridge_values[i]));
                }
            }

            for (int i = (num_divisions - 1); i > 0; i--) {
                Eigen::MatrixXd independent_vars(num_trials, 1);
                Eigen::MatrixXd dependent_vars(num_trials, 1);

                int num_paths = 0;

                // Initializing next Brownian Bridge values since they're needed for the choice of values
                for (int l = 0; l < num_trials; l++) {
                    next_brownian_bridge_values[l] = (i * T / (i + 1)) * brownian_bridge_values[l] + sigma * sqrt(i * pow(T, 2) * dt / (i + 1)) * get_gaussian();
                }

                // Choosing which values to use for regression/interpolation
                // Adapting this to non-zero cashflow. If current rebate + discounted cashflow up to exercise < discounted min rebate
                // at next exercise date, then do not evaluate. This is changed from if K - S_i > 0 then keep.
                // This is basically just another method of sub-optimal point elimination that provides a greater lower bound.

                // Taken from last sentence of first paragraph of page 12 of:
                // https://papers.ssrn.com/sol3/papers.cfm?abstract_id=1331904
                for (int j = 0; j < num_trials - 1; j++) {
                    double current_value = S_0 * exp(brownian_bridge_values[j]);
                    double next_value = S_0 * exp(next_brownian_bridge_values[j]);

                    if (call_flag == 1) {
                        if ((current_value - K) + (next_value - K) / R * exp(r * i / num_divisions * T) < 
                            (next_value - K) / R * exp(r * (i + 1) / num_divisions * T)) {
                            num_paths++;
                            independent_vars(num_paths, 0) = current_value;
                            dependent_vars(num_paths, 0) = value[j] / R;
                        }
                    } else {
                        if ((K - current_value) + (K - next_value) / R * exp(r * i / num_divisions * T) > 
                            (K - next_value) / R * exp(r * (i + 1) / num_divisions * T)) {
                            num_paths++;
                            independent_vars(num_paths, 0) = current_value;
                            dependent_vars(num_paths, 0) = value[j] / R;
                        }
                    }
                }
                
                if (num_paths > 0) {
                    // regressing the dependent_variables on the independent variables
                    int poly_degree = min(num_paths, 5);
                    Eigen::MatrixXd a_optimal(poly_degree, 1);
                    Eigen::MatrixXd X;
                    tie(X, a_optimal) = polynomial_regression(independent_vars, dependent_vars, poly_degree, num_paths, basis);

                    // Calculating the polynomial at the given point
                    for (int j = 0; j < num_trials; j++) {
                        double optimal_poly_eval = 0;
                        double current_value = S_0 * exp(brownian_bridge_values[j]);

                        // Calculating polynomial evaluation
                        if (basis == "Power") {
                            for (int l = 0; l < poly_degree; l++) {
                                // Polynomial evaluation with different bases
                                optimal_poly_eval += a_optimal(l, 0) * pow(current_value, l);
                            }
                        } else if (basis == "Laguerre") {
                            for (int l = 0; l < poly_degree; l++) {
                                double poly_eval = 0;

                                for (int m = 0; m < l + 1; m++) {
                                    poly_eval += n_choose_k(l, m) / tgamma(m + 1) * pow(-current_value, m);
                                }
                                optimal_poly_eval += poly_eval;
                            }
                        } else if (basis == "Hermitian") {
                            for (int l = 0; l < poly_degree; l++) {
                                double poly_eval = 0;

                                for (int m = 0; m < (l / 2) + 1; m++) {
                                    poly_eval += pow(-1, m) * pow(2 * current_value, l - 2 * m) / (tgamma(m + 1) * tgamma(static_cast<int>(l - 2 * m) + 1));
                                }
                                optimal_poly_eval += tgamma(l + 1) * poly_eval;
                            }
                        }


                        if (call_flag == 1) {
                            // (S - K) > poly_eval is the Andersen trigger method for which convergence is sped up.
                            // This maximizes the average cutoff over the simulated paths.
                            if (((current_value - K) > optimal_poly_eval) && ((current_value - K) > 0.0)) {
                                value[j] = current_value - K;
                            } else {
                                value[j] /= R;
                            }
                        } else {
                            if (((K - current_value) > optimal_poly_eval) && ((K - current_value) > 0.0)) {
                                value[j] = K - current_value;
                            } else {
                                value[j] /= R;
                            }
                        }
                    }

                    // Updating path values
                    for (int l = 0; l < num_trials; l++) {
                        brownian_bridge_values[l] = next_brownian_bridge_values[l];
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
	
	return 0;
}