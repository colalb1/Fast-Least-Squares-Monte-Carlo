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
#include <eigen3/Eigen/Core>

// Declaring namespaces
using Eigen::MatrixXd;
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
    for (int j = 1; j <= order; j++) {
        for (int i = 1; i <= num_obs; i++) {
            X(i, j) = pow(independent(i, 1), j - 1);
        }
    }	

	for (int i = 1; i <= num_obs; i++) {
        y(i, 1) = dependent(i, 1);
    }
	
    // solving
	Eigen::MatrixXd X_TX(order, order);
	X_TX = X.transpose() * X;

	return X_TX.inverse() * X.transpose() * y;
}

int main (int argc, char* argv[]) {
	for (int i = 0; i < 10; i++) {
        cout << get_gaussian() << "\n";
    }
    return 0;
}