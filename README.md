# Fast Least Squares Monte Carlo
Fast implementation of Least-Squares Monte Carlo aka Longstaff Schwartz in C++. This project is in development, meaning the code and description are not complete.

## Files
[fast-lsmc.cpp](https://github.com/colalb1/Fast-Least-Squares-Monte-Carlo/blob/main/fast-lsmc.cpp) contains the optimized methods for Longstaff-Schwartz and serves as a calculator.

[fast-lsmc.hpp](https://github.com/colalb1/Fast-Least-Squares-Monte-Carlo/blob/main/fast-lsmc.hpp) is the header file for [generate-data.cpp](https://github.com/colalb1/Fast-Least-Squares-Monte-Carlo/blob/main/generate-data.cpp) that contains the "function" version of the *fast-lsmc.cpp* method and a modified version of the original method provided by [R.S.](https://rsree.ise.illinois.edu/Prof._R.S._Sreenivas_%28Main%29.html) that can calculate put **and** call option prices.

[generate-data.cpp](https://github.com/colalb1/Fast-Least-Squares-Monte-Carlo/blob/main/generate-data.cpp) generates the option prices over an array of parameters and saves it as a .csv file. This has yet to be generated since I have not checked whether the method works. This needs to be optimized **badly**. These prices will be compared to the empirical price to check which is most efficient. This also tracks computation time and memory.

## Explanation

The original Longstaff Schwartz implementation was provided by [R.S.](https://rsree.ise.illinois.edu/Prof._R.S._Sreenivas_%28Main%29.html) of the UIUC Industrial Engineering Department. I will be augmenting this method.

Attempting to implement a method that chooses the best basis between Power, Laguerre, and Hermite based on the greatest $R^2_{adj}$ value via [this](https://www.sciencedirect.com/science/article/pii/S0165188913000493) paper. Adapted path conditions to non-zero cash flow to increase accuracy. This provides a greater lower bound for sub-optimal point elimination (layman: throws out more points to increase accuracy). Implemented Andersen trigger method to speed up convergence and maximize the average cutoff over simulated paths.

Might abandon choosing the best basis since the point of this is to be faster; this makes it slower although more accurate. *Edit:* The basis optimizer was not more accurate. It should have been in theory. Changed the basis to an input between the Power, Laguerre, and Hermitian bases since this decreases computation time (which was part of the point of this program anyway).

Restricted policy iteration by stratifying each potential path by its value and used continuation of previous regression.

Stratification and double-regression enhancement from [this](https://www.sciencedirect.com/science/article/pii/S0165188913000493) cannot be done since each region is only one timestep long since this is an American option as opposed to a Bermudan option. This means I must implement other optimizations. Batched iteration also prevents this since I need access to all of the paths.

Currently implementing Brownian Bridge method to reduce space requirements since only one time iteration needs to be stored at a time instead of the whole simulation. I will write a more technical explanation of this after I finish the project, but for now reference [this](https://en.wikipedia.org/wiki/Brownian_bridge) for an explanation. Simply put, I choose the last price first then walk backward toward the original price. I reckon this is not as "random" but achieves a similar outcome with much less memory since I only need to keep track of the current simulation instead of all the simulation timesteps.


## My Review of C++
Fast. Really slow to debug. This was my first time writing C++ so if there are any glaring errors or easy optimizations that I missed, let me know.
