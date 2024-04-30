// Importing packages and header file
#include <iostream>
#include <iomanip>

#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <chrono>
#include "fast-lsmc.hpp"

int main() {
    // testing if header file works
    cout << option_value(10.0, 0.1, 0.4, 100, 150, 1000, 10000, 0, "Hermitian") << endl;


    // NEED TO WRITE CSV FILE INPUT. SEE CHATGPT CSV FILE PROCESSING CODE
    return 0;
}