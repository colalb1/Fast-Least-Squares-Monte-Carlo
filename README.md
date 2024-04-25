# Fast Least Squares Monte Carlo
Fast implementation of Least-Squares Monte Carlo aka Longstaff Schwartz in C++. This project is in early development.

The original Longstaff Schwartz implementation was provided by [R.S.](https://rsree.ise.illinois.edu/Prof._R.S._Sreenivas_%28Main%29.html) of the UIUC Industrial Engineering Department. I will be augmenting this method.

Attempting to implement a method that chooses the best basis between Power, Laguerre, and Hermite based on the greatest $R^2_{adj}$ value via [this](https://www.sciencedirect.com/science/article/pii/S0165188913000493) paper. Adapted path conditions to non-zero cash flow to increase accuracy. This provides a greater lower bound for sub-optimal point elimination (layman: throws out more points to increase accuracy). Implemented Andersen trigger method to speed up convergence and maximize the average cutoff over simulated paths.

Might abandon choosing the best basis since the point of this is to be faster; this makes it slower although more accurate.


## My Review of C++
Fast. Also makes me want to scoop my eyes out.
