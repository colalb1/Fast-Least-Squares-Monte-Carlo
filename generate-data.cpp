// Importing packages and header file
#include <iostream>
#include <iomanip>

#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <chrono>
#include "fast-lsmc.hpp"

using namespace std;

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
    vector<double> basis_store = {};
    vector<double> call_flag_store = {};

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
                    double temp_price = option_value_original(expiry_time[i], 
                                                              risk_free_rates[j], 
                                                              volatilities[k], 
                                                              starting_price, 
                                                              call_strikes[l], 
                                                              num_divisions_, 
                                                              num_sims_, 
                                                              1);

                    // Bases
                    for (int m = 0; m < sizeof(bases); m++) {
                        // Saving data
                        expiry_time_store.push_back(expiry_time[i]);
                        risk_free_rates_store.push_back(risk_free_rates[j]);
                        volatilities_store.push_back(volatilities[k]);
                        strikes_store.push_back(call_strikes[l]);
                        call_flag_store.push_back(1);
                        basis_store.push_back(m);

                        original_price_store.push_back(temp_price);
                        optimized_price_store.push_back(option_value(expiry_time[i],
                                                                     risk_free_rates[j], 
                                                                     volatilities[k], 
                                                                     starting_price, 
                                                                     call_strikes[l], 
                                                                     num_divisions_, 
                                                                     num_sims_, 
                                                                     1, 
                                                                     bases[m]));
                    }
                }
                // Puts
                for (int l = 0; l < put_strikes.size(); l++) {
                    // Price of the original method since I will be saving it for each basis function
                    double temp_price = option_value_original(expiry_time[i], 
                                                              risk_free_rates[j], 
                                                              volatilities[k], 
                                                              starting_price, 
                                                              put_strikes[l], 
                                                              num_divisions_, 
                                                              num_sims_, 
                                                              0);

                    // Bases
                    for (int m = 0; m < sizeof(bases); m++) {
                        // Saving data
                        expiry_time_store.push_back(expiry_time[i]);
                        risk_free_rates_store.push_back(risk_free_rates[j]);
                        volatilities_store.push_back(volatilities[k]);
                        strikes_store.push_back(put_strikes[l]);
                        call_flag_store.push_back(0);
                        basis_store.push_back(m);

                        original_price_store.push_back(temp_price);
                        optimized_price_store.push_back(option_value(expiry_time[i],
                                                                     risk_free_rates[j], 
                                                                     volatilities[k], 
                                                                     starting_price, 
                                                                     put_strikes[l], 
                                                                     num_divisions_, 
                                                                     num_sims_, 
                                                                     0, 
                                                                     bases[m]));
                    }
                } 
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

    write_csv("option_pricing.csv", output);
    return 0;
}