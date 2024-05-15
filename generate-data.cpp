// Importing packages and header file
#include <iostream>
#include <iomanip>

#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <chrono>
#include "fast-lsmc.hpp"
#include <thread>

using namespace std;

// Standard normal cdf
double standard_normal_cdf(double x) {
    return (1 + erf(x / sqrt(2))) / 2;
}

// Black-Scholes empirical option price
// Using this for error comparison
double empirical_option(double S_0, double T, double sigma, double K, double r, int call_flag) {
    double d_1 = (log(S_0 / K) + (r + pow(sigma, 2) / 2) * T) / (sigma * sqrt(T));
    double d_2 = d_1 - sigma * sqrt(T);

    if (call_flag == 1) {
        return S_0 * standard_normal_cdf(d_1) - K * exp(-r * T) * standard_normal_cdf(d_2);
    }
    return K * exp(-r * T) * standard_normal_cdf(-d_2) - S_0 * standard_normal_cdf(-d_1);
}

// The function below is a modified version of this: https://www.gormanalysis.com/blog/reading-and-writing-csv-files-with-cpp/
void write_csv(string filename, vector<pair<string, vector<double>>> dataset) {
    // Make a CSV file with one or more columns of double values
    // Each column of data is represented by the pair <column name, column data>
    //   as pair<string, vector<double>>
    // The dataset is represented as a vector of these columns
    // Note that all columns should be the same size
    
    // Create an output filestream object
    ofstream myFile(filename);
    
    // Send column names to the stream
    for (int j = 0; j < dataset.size(); ++j) {
        myFile << dataset.at(j).first;
        if (j != dataset.size() - 1) {
            myFile << ","; // No comma at end of line
        }
    }
    myFile << "\n";
    
    // Send data to the stream
    for (int i = 0; i < dataset.at(0).second.size(); ++i) {
        for (int j = 0; j < dataset.size(); ++j) {
            myFile << dataset.at(j).second.at(i);
            if (j != dataset.size() - 1) {
                myFile << ","; // No comma at end of line
            }
        }
        myFile << "\n";
    }
    // Close the file
    myFile.close();
}


int main() {
    // Using num_divisions = 2000, num_sims = 10000, and S_0 = 100 for all.
    // I would check more values but I am "computationally limited". My computer is 5 years old.

    // Initializing parameters
    vector<double> expiry_time = {1, 5, 10};
    vector<double> risk_free_rates = {0.02, 0.03, 0.05, 0.1};
    vector<double> volatilities = {0.1, 0.2, 0.3, 0.4, 0.5};
    int starting_price = 100;
    int num_divisions_ = 2000;
    int num_sims_ = 10000;
    string bases[3] = {"Power", "Hermitian", "Laguerre"}; // {0, 1, 2} for encoding purposes
    vector<double> call_strikes = {102, 105, 110, 120, 150};
    vector<double> put_strikes = {98, 95, 90, 80, 50};

    // Initializing storage for data table
    vector<double> expiry_time_store = {};
    vector<double> risk_free_rates_store = {};
    vector<double> volatilities_store = {};
    vector<double> strikes_store = {};
    vector<double> original_price_store = {};
    vector<double> optimized_price_store = {};
    vector<double> empirical_price = {};
    vector<double> basis_store = {};
    vector<double> call_flag_store = {};
    vector<double> original_computation_time = {};
    vector<double> optimized_computation_time = {};
    

    // creating mutex to protect the vectors I push back on
    mutex vector_mutex;


    // Excuse the asymptotic runtime of this for loop. This is a time for which Python would exponentially simplify things.
    // Calculating prices and filling in data vectors

    // Expiry times
    for (int i = 0; i < expiry_time.size(); i++) {
        // Risk-free rates
        for (int j = 0; j < risk_free_rates.size(); j++) {
            // Volatilities
            for (int k = 0; k < volatilities.size(); k++) {
                // Calls
                for (int l = 0; l < call_strikes.size(); l++) {
                    // Price of the original method since I will be saving it for each basis function
                    // Starting computation time
                    const clock_t original_begin_time = clock();
                    double temp_price = option_value_original(expiry_time[i], 
                                                              risk_free_rates[j], 
                                                              volatilities[k], 
                                                              starting_price, 
                                                              call_strikes[l], 
                                                              num_divisions_, 
                                                              num_sims_, 
                                                              1);
                    // Stopping computation time and calculating
                    double original_duration = double(clock() - original_begin_time) / CLOCKS_PER_SEC;

                    // Empirical price according to Black-Scholes
                    double temp_empirical_price = empirical_option(starting_price, 
                                                                   expiry_time[i], 
                                                                   volatilities[k], 
                                                                   call_strikes[l], 
                                                                   risk_free_rates[j],
                                                                   1);
                    

                    // Setting up multithreading thread ppol
                    vector<thread> thread_pool;


                    // lambda functions for each basis computation. This is verbose and NEEDS to be optimized.
                    auto compute_task = [&expiry_time_store, &expiry_time, &risk_free_rates_store, &risk_free_rates,
                                         &volatilities_store, &volatilities, &strikes_store, &call_strikes,
                                         &call_flag_store, &basis_store, &original_price_store, &original_computation_time,
                                         &optimized_price_store, &starting_price, &num_divisions_, &num_sims_, &bases, 
                                         &i, &j, &k, &l, &temp_price, &original_duration, 
                                         &optimized_computation_time, &vector_mutex, 
                                         &empirical_price, &temp_empirical_price](int basis_index) {
                        // Synchronization mutex lock
                        lock_guard<mutex> lock(vector_mutex);

                        // Saving data
                        expiry_time_store.push_back(expiry_time[i]);
                        risk_free_rates_store.push_back(risk_free_rates[j]);
                        volatilities_store.push_back(volatilities[k]);
                        strikes_store.push_back(call_strikes[l]);
                        call_flag_store.push_back(1);
                        basis_store.push_back(basis_index);
                        empirical_price.push_back(temp_empirical_price);

                        original_price_store.push_back(temp_price);
                        original_computation_time.push_back(original_duration);

                        const clock_t optimized_begin_time = clock();
                        optimized_price_store.push_back(option_value(expiry_time[i],
                                                                     risk_free_rates[j], 
                                                                     volatilities[k], 
                                                                     starting_price, 
                                                                     call_strikes[l], 
                                                                     num_divisions_, 
                                                                     num_sims_, 
                                                                     1, 
                                                                     bases[basis_index]));
                        double optimized_duration = double (clock() - optimized_begin_time) / CLOCKS_PER_SEC;
                        optimized_computation_time.push_back(optimized_duration);
                    };

                    // Enqueue computations to thread pool
                    for (int m = 0; m < 3; m++) {           // 3 = sizeof(bases)
                        lock_guard<mutex> lock(vector_mutex);

                        // Want to save the function call, but cannot call it straight into the thread due to the
                        // interaction between the two
                        function<void()> temp_compute = bind(compute_task, m);
                        thread_pool.emplace_back(thread(temp_compute));
                    }

                    // Join threads to main thread
                    for (auto& thread:thread_pool) {
                        if (thread.joinable()) {
                            thread.join();
                        }
                    }
                }
                
                // Puts
                for (int l = 0; l < put_strikes.size(); l++) {
                    // Price of the original method since I will be saving it for each basis function
                    // Starting computation time
                    const clock_t original_begin_time = clock();
                    double temp_price = option_value_original(expiry_time[i], 
                                                              risk_free_rates[j], 
                                                              volatilities[k], 
                                                              starting_price, 
                                                              put_strikes[l], 
                                                              num_divisions_, 
                                                              num_sims_, 
                                                              0);
                    // Stopping computation time and calculating
                    double original_duration = double(clock() - original_begin_time) / CLOCKS_PER_SEC;

                    // Empirical price according to Black-Scholes
                    double temp_empirical_price = empirical_option(starting_price, 
                                                                   expiry_time[i], 
                                                                   volatilities[k], 
                                                                   call_strikes[l], 
                                                                   risk_free_rates[j],
                                                                   0);

                    // Setting up multithreading thread ppol
                    vector<thread> thread_pool;


                    // lambda functions for each basis computation. This is verbose and NEEDS to be optimized.
                    auto compute_task = [&expiry_time_store, &expiry_time, &risk_free_rates_store, &risk_free_rates,
                                         &volatilities_store, &volatilities, &strikes_store, &put_strikes,
                                         &call_flag_store, &basis_store, &original_price_store, &original_computation_time,
                                         &optimized_price_store, &starting_price, &num_divisions_, &num_sims_, &bases, 
                                         &i, &j, &k, &l, &temp_price, &original_duration, 
                                         &optimized_computation_time, &vector_mutex,
                                         &empirical_price, &temp_empirical_price](int basis_index) {
                        // Synchronization mutex lock
                        lock_guard<mutex> lock(vector_mutex);

                        // Saving data
                        expiry_time_store.push_back(expiry_time[i]);
                        risk_free_rates_store.push_back(risk_free_rates[j]);
                        volatilities_store.push_back(volatilities[k]);
                        strikes_store.push_back(put_strikes[l]);
                        call_flag_store.push_back(0);
                        basis_store.push_back(basis_index);
                        empirical_price.push_back(temp_empirical_price);

                        original_price_store.push_back(temp_price);
                        original_computation_time.push_back(original_duration);

                        const clock_t optimized_begin_time = clock();
                        optimized_price_store.push_back(option_value(expiry_time[i],
                                                                     risk_free_rates[j], 
                                                                     volatilities[k], 
                                                                     starting_price, 
                                                                     put_strikes[l], 
                                                                     num_divisions_, 
                                                                     num_sims_, 
                                                                     0, 
                                                                     bases[basis_index]));
                        double optimized_duration = double(clock() - optimized_begin_time) / CLOCKS_PER_SEC;
                        optimized_computation_time.push_back(optimized_duration);
                    };

                    // Enqueue computations to thread pool
                    for (int m = 0; m < 3; m++) {           // 3 = sizeof(bases)
                        lock_guard<mutex> lock(vector_mutex);

                        // Want to save the function call, but cannot call it straight into the thread due to the
                        // interaction between the two
                        function<void()> temp_compute = bind(compute_task, m);
                        thread_pool.emplace_back(thread(temp_compute));
                    }

                    // Join threads to main thread
                    for (auto& thread:thread_pool) {
                        if (thread.joinable()) {
                            thread.join();
                        }
                    }
                }
                cout << "Finished:" << endl;
                cout << volatilities[k] << endl;
                cout << risk_free_rates[j] << endl;
                cout << expiry_time[i] << endl;
                cout << endl;
            }
        }
    }
    // Forming the data table
    vector<pair<string, vector<double>>> output;
    
    output.push_back({"Expiry Time", expiry_time_store});
    output.push_back({"Risk-Free Rate", risk_free_rates_store});
    output.push_back({"Volatility", volatilities_store});
    output.push_back({"Strike Price", strikes_store});
    output.push_back({"Original Method Price", original_price_store});
    output.push_back({"Optimized Method Price", optimized_price_store});
    output.push_back({"Basis", basis_store});
    output.push_back({"Call Flag", call_flag_store});
    output.push_back({"Original Computation Time", original_computation_time});
    output.push_back({"Optimized Computation Time", optimized_computation_time});

    write_csv("option_pricing.csv", output);
    return 0;
}