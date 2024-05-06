# Fast Least Squares Monte Carlo
Fast implementation of Least-Squares Monte Carlo (Longstaff-Schwartz) for American options in C++. This project is in development; the code and description are incomplete.

## Precursor

The original Longstaff-Schwartz implementation was provided by [R.S.](https://rsree.ise.illinois.edu/Prof._R.S._Sreenivas_%28Main%29.html) of the UIUC Industrial Engineering Department. The genesis of augmenting this method came from discussing option pricing using Monte Carlo methods after class.

## Files
[*fast-lsmc.cpp*](https://github.com/colalb1/Fast-Least-Squares-Monte-Carlo/blob/main/fast-lsmc.cpp) contains the optimized methods for Longstaff-Schwartz and serves as a calculator.

[*fast-lsmc.hpp*](https://github.com/colalb1/Fast-Least-Squares-Monte-Carlo/blob/main/fast-lsmc.hpp) is the header file for [*generate-data.cpp*](https://github.com/colalb1/Fast-Least-Squares-Monte-Carlo/blob/main/generate-data.cpp) that contains the "function" version of the *fast-lsmc.cpp* method and a modified version of the original method provided by [R.S.](https://rsree.ise.illinois.edu/Prof._R.S._Sreenivas_%28Main%29.html) that can calculate put **and** call option prices.

[*generate-data.cpp*](https://github.com/colalb1/Fast-Least-Squares-Monte-Carlo/blob/main/generate-data.cpp) generates the option prices over an array of parameters and saves it as a .csv file. The estimated prices and computation times will be used to compare the efficiency of each method.

## Optimizations

These optimizations are meant to make the original Longstaff-Schwartz more efficient. The following contains a short explanation of each improvement.

### Choosing the Basis

After reviewing literature, choosing a basis other than the standard power basis (that was used in the [original Longstaff-Schwartz paper](https://people.math.ethz.ch/~hjfurrer/teaching/LongstaffSchwartzAmericanOptionsLeastSquareMonteCarlo.pdf)) was consistent amongst many papers. The three tested bases were the Power, Hermitian, and Laguerre bases. They are defined as follows (**k** is the number of desired basis functions):

**Power:** $$\sum_{n = 0}^k x^n$$

**Hermitian:** $$\sum_{n = 0}^k n!\sum_{m = 0}^{\lfloor \frac{n}{2} \rfloor} (-1)^{m}(2x)^{n - 2m} * \frac{1}{m!(n - 2m)!}$$

**Laguerre:** $$\sum_{n = 0}^k \sum_{m = 0}^n \frac{(-x)^m}{m!} {n \choose m}$$

I will omit further explanation of the bases for brevity; refer to the top of page 6 of [this](https://jfin-swufe.springeropen.com/articles/10.1186/s40854-015-0019-0) paper for more clarity regarding basis construction.

I implemented a method that programmatically chose the optimal basis based on which had the greatest $R^2_{adj}$ value via [this](https://www.sciencedirect.com/science/article/pii/S0165188913000493) paper. I abandoned this idea because there is no way (that I know of) around computing the regression three times for each iteration, and this significantly increased computation time. This works well in theory, but having nearly three times the computation time for (even) a significant accuracy gain is not worth it.

### Path Conditions

In Least-Squares Monte Carlo methods for options pricing, the condition $(S - K) > 0$ or $(K - S) > 0$ for calls and puts, respectively, is imposed on the price paths to select the path data that will be used for regression to generate the equation that predicts the option price at the previous step. [This](https://www.sciencedirect.com/science/article/pii/S0165188913000493) paper tightened this restriction to a non-zero cash flow and showed an accuracy increase due to the greater lower bound for sub-optimal point eliminations. Thus, the condition changes to the following (assume time $t$ is discrete):

$$(S_{t} - K) + (S_{t + 1} - K) * \exp(-r(T - 2t)) < (S_{t + 1} - K) * \exp(-r(T - 2(t + 1)))$$

for calls and

$$(K - S_{t}) + (K - S_{t + 1}) * \exp(-r(T - 2t)) > (K - S_{t + 1}) * \exp(-r(T - 2(t + 1)))$$

for puts. Put simply, a path is excluded if the price in the previous step does not reach a certain threshold based on the decay imposed by the risk-free rate.

The Andersen trigger method (otherwise known as LSA) was also added to improve the convergence rate by maximizing the average cutoff over simulated paths. This is essentially adding a constant to one side of the above inequality. More details can be found on page 12 of the paper linked above.


### Policy Iteration

Restricted policy iteration by stratifying each potential path by its value and used continuation of previous regression.

### Stratification and Double-Regression Enhancement

Stratification and double-regression enhancement from [this](https://www.sciencedirect.com/science/article/pii/S0165188913000493) cannot be done since each region is only one timestep long since this is an American option as opposed to a Bermudan option. This means I must implement other optimizations. Batched iteration also prevents this since I need access to all of the paths.


### Brownian Bridge 

Currently implementing Brownian Bridge method to reduce space requirements since only one time iteration needs to be stored at a time instead of the whole simulation. I will write a more technical explanation of this after I finish the project, but for now reference [this](https://en.wikipedia.org/wiki/Brownian_bridge) for an explanation. Simply put, I choose the last price first then walk backward toward the original price. I reckon this is not as "random" but achieves a similar outcome with much less memory since I only need to keep track of the current simulation instead of all the simulation timesteps.

## Conclusion

Add metrics on performance and discuss improvements in general.

## My Review of C++
Fast. Really slow to debug. This was my first time writing C++ so if there are any glaring errors or easy optimizations that I missed, let me know.
