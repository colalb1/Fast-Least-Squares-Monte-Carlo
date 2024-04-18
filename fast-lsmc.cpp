#include <iostream>
#include <iomanip>

#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <random>
#include <chrono>

using namespace std;

unsigned seed = (unsigned) std::chrono::system_clock::now().time_since_epoch().count();
default_random_engine generator;

double risk_free_rate, strike_price, initial_stock_price, expiration_time, volatility;
int no_of_trials, no_of_divisions, no_of_simulations;

double max(double a, double b) {
	return (b < a) ? a:b;
}

int min(int a, int b) {
	return (b < a) ? b:a;
}

double get_uniform() {
    // http://www.cplusplus.com/reference/random/exponential_distribution/
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    double number = distribution(generator);
    return number;
}

double get_gaussian() {
	return sqrt(-2.0 * log(get_uniform())) * cos(2 * M_PI * get_uniform());
}

int main (int argc, char* argv[]) {
	for (int i = 0; i < 10; i++) {
        cout << i << "\n";
    }
    return 0;
}
