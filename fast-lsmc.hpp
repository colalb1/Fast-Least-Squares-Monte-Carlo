#ifndef MY_HEADER_H
#define MY_HEADER_H

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
int num_trials, num_divisions, num_sims, call_flag, no_of_trials;
string basis;

// Max and min functions
static double max(double a, double b) {
	return (b < a) ? a:b;
}

static int min(int a, int b) {
	return (b < a) ? b:a;
}

// Generating standard uniform
static double get_uniform() {
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    return distribution(generator);
}

// Box-Muller transform
static double get_gaussian() {
	return sqrt(-2.0 * log(get_uniform())) * cos(2 * M_PI * get_uniform());
}

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
	Eigen::MatrixXd X(num_obs, order);
	Eigen::MatrixXd y(num_obs, 1);
    
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

double option_value(double T, double r, double sigma, double S_0, double K, int num_divisions, int num_sims, int call_flag, string basis) {
    // defining time, rate, and SD step lengths
    double dt = T / ((double) num_divisions);
	double dR = (r - 0.5 * pow(sigma, 2)) * dt;
	double dSD = sigma * sqrt(dt);
	double R = exp(r * T / ((double) num_divisions));
    
    // defining American put price
    double option_price = 0.0;
    
    for (int k = 0; k <= (num_sims / 200); k++) {
		
		// setting trial number
		if (k != num_sims / 200) {
            num_trials = 200;
        } else {
            num_trials = num_sims % 200;
        }			
		
		if (num_trials > 0) {
			double brownian_bridge_values[num_trials], next_brownian_bridge_values[num_trials], value[num_trials];
            for (int i = 0; i < num_trials; i++) {
                brownian_bridge_values[i] = (r - pow(sigma, 2) / 2) * T + sigma * sqrt(T) * get_gaussian();

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
                    int poly_degree = min(num_paths, 5);
                    Eigen::MatrixXd a_optimal(poly_degree, 1);
                    Eigen::MatrixXd X;
                    tie(X, a_optimal) = polynomial_regression(independent_vars, dependent_vars, poly_degree, num_paths, basis);

                    for (int j = 0; j < num_trials; j++) {
                        double optimal_poly_eval = 0;
                        double current_value = S_0 * exp(brownian_bridge_values[j]);

                        if (basis == "Power") {
                            for (int l = 0; l < poly_degree; l++) {
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
	
	if (num_sims % 200 == 0) {
        return option_price / ((double) num_sims / 200);
    }
	
    return option_price / ((double) num_sims / 200 + 1);
}


// Code below this is not originally mine (modified slightly to generalize and fit to this application)


Eigen::MatrixXd polynomial_regression_original(Eigen::MatrixXd Independent_Variables, Eigen::MatrixXd Dependent_Variable, int order, int no_of_observations)
{
	Eigen::MatrixXd X(no_of_observations, order);
	Eigen::MatrixXd Y(no_of_observations, 1);
	
	for (int i = 1; i <= no_of_observations; i++)
		Y(i,1) = Dependent_Variable(i,1);
	
	for (int j = 1; j <= order; j++) 
		for (int i = 1; i <= no_of_observations; i++) 
			X(i,j) = pow(Independent_Variables(i,1), j-1);
	
	// return inv(XT*X)*XT*Y
	Eigen::MatrixXd X_transpose_times_X(order, order);
	X_transpose_times_X = X.transpose()*X;
	return (X_transpose_times_X.inverse() * X.transpose() * Y);
}

double option_value_original(double expiration_time, double risk_free_rate, double volatility, double initial_stock_price, double strike_price, int no_of_divisions, int no_of_simulations, int call_flag) {
	double delta_T = expiration_time/((double) no_of_divisions);
	double delta_R = (risk_free_rate - 0.5*pow(volatility,2))*delta_T;
	double delta_SD = volatility*sqrt(delta_T);
	double R = exp(risk_free_rate*expiration_time/((double) no_of_divisions));

	
	// given array size limitations, I will run batches of 200 runs if no_of_simulations 
	// exceeds 200.  
	
	double option_price = 0.0;
	for (int k = 0; k < (no_of_simulations/200 + 1); k++) {
		
		if (k != no_of_simulations/200) 
			no_of_trials = 200;
		else 
			no_of_trials = no_of_simulations%200;
		
		if (no_of_trials != 0) {
			double asset_price[200][no_of_divisions];
			
			for (int i = 0; i < no_of_trials; i++)
				asset_price[i][0] = initial_stock_price;
			
			for (int i = 0; i < no_of_trials; i++) 
				for (int j = 1; j < no_of_divisions; j++) 
					asset_price[i][j] = asset_price[i][j-1]*exp(delta_R + delta_SD*get_gaussian());

			
			double value[no_of_trials];
			// initialize the value based on the price at final stage
			for (int i = 0; i < no_of_trials; i++) 
				if (call_flag == 0)
					value[i] = max(0.0, strike_price - asset_price[i][no_of_divisions-1]);
				else
					value[i] = max(0.0, asset_price[i][no_of_divisions-1] - strike_price);
			
			for (int i = (no_of_divisions-1); i > 0; i--) {
				Eigen::MatrixXd independent_variables(no_of_trials,1);
				Eigen::MatrixXd dependent_variables(no_of_trials,1);
				int no_of_variables = 0;
				for (int j = 0; j < no_of_trials; j++) {
					if (max(0.0, strike_price - asset_price[j][i]) > 0 && call_flag == 0) {
						no_of_variables++;
						independent_variables(no_of_variables, 0) = asset_price[j][i];
						dependent_variables(no_of_variables, 0) = value[j]/R;
					} else if (max(0.0, asset_price[j][i] - strike_price) > 0 && call_flag == 1) {
						no_of_variables++;
						independent_variables(no_of_variables, 0) = asset_price[j][i];
						dependent_variables(no_of_variables, 0) = value[j]/R;
					}
				}
                
				
				if (no_of_variables > 0) {
					// regressing the dependent_variables on the independent variables using a 4th order polynomial
					Eigen::MatrixXd a(min(5,no_of_variables),1);
					a = polynomial_regression_original(independent_variables, dependent_variables, min(5,no_of_variables), no_of_variables);
					if (no_of_variables >= 5) {
						for (int j = 0; j < no_of_trials; j++) {
							if (call_flag == 0) {
                                if ( ((strike_price - asset_price[j][i]) > (a(0,0) + 
																			(a(1,0)*asset_price[j][i]) + 
																			(a(2,0)*pow(asset_price[j][i],2)) + 
																			(a(3,0)*pow(asset_price[j][i],3)) + 
																			(a(4,0)*pow(asset_price[j][i],4)))) && 
									( (strike_price -asset_price[j][i]) > 0.0 ) )
									value[j] = strike_price - asset_price[j][i];
								else
									value[j] = value[j]/R;
                            } else {
                                if ( ((asset_price[j][i] - strike_price) > (a(0,0) + 
																			(a(1,0)*asset_price[j][i]) + 
																			(a(2,0)*pow(asset_price[j][i],2)) + 
																			(a(3,0)*pow(asset_price[j][i],3)) + 
																			(a(4,0)*pow(asset_price[j][i],4)))) && 
									( (asset_price[j][i] - strike_price) > 0.0 ) )
									value[j] = asset_price[j][i] - strike_price;
								else
									value[j] = value[j]/R;
                            }
						}
					}
					else if (no_of_variables == 4) {
						for (int j = 0; j < no_of_trials; j++) {
							if (call_flag == 0)
								if ( ((strike_price - asset_price[j][i]) > (a(0,0) + 
																			(a(1,0)*asset_price[j][i]) + 
																			(a(2,0)*pow(asset_price[j][i],2)) + 
																			(a(3,0)*pow(asset_price[j][i],3)))) && 
									( (strike_price -asset_price[j][i]) > 0.0 ) )
									value[j] = strike_price - asset_price[j][i];
								else
									value[j] = value[j]/R;
							else
								if ( ((asset_price[j][i] - strike_price) > (a(0,0) + 
																			(a(1,0)*asset_price[j][i]) + 
																			(a(2,0)*pow(asset_price[j][i],2)) + 
																			(a(3,0)*pow(asset_price[j][i],3)))) && 
									( (asset_price[j][i] - strike_price) > 0.0 ) )
									value[j] = asset_price[j][i] - strike_price;
								else
									value[j] = value[j]/R;
						}
						
					}
					else if (no_of_variables == 3) {
						for (int j = 0; j < no_of_trials; j++) {
							if (call_flag == 0)
								if ( ((strike_price - asset_price[j][i]) > (a(0,0) + 
																			(a(1,0)*asset_price[j][i]) + 
																			(a(2,0)*pow(asset_price[j][i],2)))) &&
									( (strike_price -asset_price[j][i]) > 0.0 ) )
									value[j] = strike_price - asset_price[j][i];
								else
									value[j] = value[j]/R;
							else
								if ( ((asset_price[j][i] - strike_price) > (a(0,0) + 
																			(a(1,0)*asset_price[j][i]) + 
																			(a(2,0)*pow(asset_price[j][i],2)))) &&
									( (asset_price[j][i] - strike_price) > 0.0 ) )
									value[j] = asset_price[j][i] - strike_price;
								else
									value[j] = value[j]/R;
						}	
					}
					else if (no_of_variables == 2) {
						for (int j = 0; j < no_of_trials; j++) {
							if (call_flag == 0)
								if ( ((strike_price - asset_price[j][i]) > (a(0,0) + 
																			(a(1,0)*asset_price[j][i]))) &&
									( (strike_price -asset_price[j][i]) > 0.0 ) )
									value[j] = strike_price - asset_price[j][i];
								else
									value[j] = value[j]/R;
							else
								if ( ((asset_price[j][i] - strike_price) > (a(0,0) + 
																			(a(1,0)*asset_price[j][i]))) &&
									( (asset_price[j][i] - strike_price) > 0.0 ) )
									value[j] = asset_price[j][i] - strike_price;
								else
									value[j] = value[j]/R;
						}
					}
					else  {
						for (int j = 0; j < no_of_trials; j++) {
							if (call_flag == 0)
								if ( ((strike_price - asset_price[j][i]) > a(0,0))  &&
									( (strike_price -asset_price[j][i]) > 0.0 ) )
									value[j] = strike_price - asset_price[j][i];
								else
									value[j] = value[j]/R;
							else
								if ( ((asset_price[j][i] - strike_price) > a(0,0))  &&
									( (asset_price[j][i] - strike_price) > 0.0 ) )
									value[j] = asset_price[j][i] - strike_price;
								else
									value[j] = value[j]/R;
						}	
					}
					
				}
			}
			
			double local_option_price = 0.0;
			for (int j = 0; j < no_of_trials; j++) 
				local_option_price += value[j];
			local_option_price = (local_option_price/((float) no_of_trials))/R;
			option_price += local_option_price;
		}
	}
	
	if (no_of_simulations%200 == 0) 
		return option_price/((double) no_of_simulations/200);
	return option_price/((double) no_of_simulations/200 + 1);
}

#endif