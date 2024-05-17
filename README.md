# Optimized Least Squares Monte Carlo
Fast implementation of Least-Squares Monte Carlo (Longstaff-Schwartz) for American options in C++. This project is in development; the code and description are incomplete. This project aims to optimize the original Longstaff-Schwartz method to greatly increase accuracy and decrease computation time given equivalent computational parameters (number of grid divisions) via mathematical and programmatic improvements.

## Precursor

The original Longstaff-Schwartz implementation was provided by [R.S.](https://rsree.ise.illinois.edu/Prof._R.S._Sreenivas_%28Main%29.html) of the UIUC Industrial Engineering Department. The idea to augment this method came from discussing option pricing using Monte Carlo methods after class.

## Files
[*fast-lsmc.cpp*](https://github.com/colalb1/Fast-Least-Squares-Monte-Carlo/blob/main/fast-lsmc.cpp) contains the optimized methods for Longstaff-Schwartz and serves as a calculator.

[*fast-lsmc.hpp*](https://github.com/colalb1/Fast-Least-Squares-Monte-Carlo/blob/main/fast-lsmc.hpp) is the header file for [*generate-data.cpp*](https://github.com/colalb1/Fast-Least-Squares-Monte-Carlo/blob/main/generate-data.cpp) that contains the "function" version of the *fast-lsmc.cpp* method and a modified version of the original method provided by [R.S.](https://rsree.ise.illinois.edu/Prof._R.S._Sreenivas_%28Main%29.html) that can calculate put **and** call option prices.

[*generate-data.cpp*](https://github.com/colalb1/Fast-Least-Squares-Monte-Carlo/blob/main/generate-data.cpp) generates the option prices over an array of parameters and saves it as a .csv file. The estimated prices and computation times will be used to compare the efficiency of each method.

[*option_pricing.csv*](https://github.com/colalb1/Fast-Least-Squares-Monte-Carlo/blob/main/option_pricing.csv) contains the data generated from the file listed above. This will be used for data analysis regarding accuracy and computation time to draw general conclusions and those regarding the bases.

[*optimization_analysis.ipynb*](https://github.com/colalb1/Fast-Least-Squares-Monte-Carlo/blob/main/optimization_analysis.ipynb) contains data analysis comparing the accuracy of the original and option pricing methods.

## Optimizations

These optimizations are meant to make the original Longstaff-Schwartz more efficient. The following contains a short explanation of each improvement.

The [eigen 3.4.0](https://eigen.tuxfamily.org/index.php?title=Main_Page) C++ package was used for all linear algebra needs.

### Choosing the Basis

After reviewing literature, choosing a basis other than the standard power basis (that was used in the [original Longstaff-Schwartz paper](https://people.math.ethz.ch/~hjfurrer/teaching/LongstaffSchwartzAmericanOptionsLeastSquareMonteCarlo.pdf)) was consistent amongst many papers. The three tested bases were the Power, Hermitian, and Laguerre bases. They are defined as follows (**k** is the number of desired basis functions):

**Power:** 

```math
\left\{x^n\right\}_{n = 0}^{k}
```

**Hermitian:**

```math
\left\{n!\sum_{m = 0}^{\lfloor \frac{n}{2} \rfloor} (-1)^{m}(2x)^{n - 2m} * \frac{1}{m!(n - 2m)!}\right\}_{n = 0}^k
```

**Laguerre:**

```math
\left\{\sum_{m = 0}^n \frac{(-x)^m}{m!} {n \choose m}\right\}_{n = 0}^k
```

I will omit further explanation of the bases for brevity; refer to the top of page 6 of [this](https://jfin-swufe.springeropen.com/articles/10.1186/s40854-015-0019-0) paper for more clarity regarding basis construction.

I implemented a method that programmatically chose the optimal basis based on which had the greatest $R^2_{adj}$ value via [this](https://www.sciencedirect.com/science/article/pii/S0165188913000493) paper. I abandoned this idea because there is no way (that I know of) around computing the regression three times for each iteration, and this significantly increased computation time. This works well in theory, but having nearly three times the computation time for a small accuracy gain is not worth it.

### Path Conditions

In Least-Squares Monte Carlo methods for options pricing, the condition $(S - K) > 0$ or $(K - S) > 0$ for calls and puts, respectively, is imposed on the price paths to select the path data that will be used for regression to generate the equation that predicts the option price at the previous step. [This](https://www.sciencedirect.com/science/article/pii/S0165188913000493) paper tightened this restriction to a non-zero cash flow and showed an accuracy increase due to the greater lower bound for sub-optimal point eliminations. Thus, the condition changes to the following (assume time $t$ is discrete):

$$(S_{t} - K) + (S_{t + 1} - K) * \exp(-r(T - 2t)) < (S_{t + 1} - K) * \exp(-r(T - 2(t + 1)))$$

for calls and

$$(K - S_{t}) + (K - S_{t + 1}) * \exp(-r(T - 2t)) > (K - S_{t + 1}) * \exp(-r(T - 2(t + 1)))$$

for puts. Put simply, a path is excluded if the price in the previous step does not reach a certain threshold based on the decay imposed by the risk-free rate.

The Andersen trigger method (otherwise known as LSA) was also added to improve the convergence rate by maximizing the average cutoff over simulated paths. This is essentially adding a constant to one side of the above inequality. More details can be found on page 12 of the paper linked above.

### Brownian Bridge 

A Brownian Bridge was simulated instead of a Brownian Motion as a Brownian Bridge requires only the last iteration of movement to be stored whereas Brownian Motion requires storage of the whole walk. This saves memory and decreases computation time since a greater portion of the cached information resides in the CPU. A Brownian Bridge is essentially when the final value of a walk is chosen first, then a pseudo-random iteration process walks the process back to some starting value (basically a Brownian Motion in reverse). I will now present the mathematical formulation of this concept.

Suppose $Z\sim N(0, 1)$ is a random sample from a standard normal distribution. Initialize the final value of the walk such that 

$$S(T) = S(0) * \exp(X(T))$$

where $S(0)$ is the starting price, and $X(T) = (r - \frac{\sigma ^ 2}{2}) + \sigma\sqrt{T} * Z$, where $\sqrt{T}Z = W(T)\sim N(0, T)$.

For each previous timestep, walk backward using the following equation:

$$X(t_i) = \frac{t_{i + 1}}{t_i}X(t_{i + 1}) + \sigma\sqrt{\frac{t_i}{t_{i + 1}} * \Delta t} * Z$$

One may observe that this achieves the starting value at $t = 0$. Therefore, the construction of the Brownian Bridge is complete.

This formulation was taken from [this](http://www.diva-portal.org/smash/get/diva2:818128/FULLTEXT01.pdf) article.


### Stratification and Double-Regression Enhancement

Another attempted optimization was the stratification and a double-regression enhancement (from [this](https://www.sciencedirect.com/science/article/pii/S0165188913000493) paper). This achieves essentially the same outcome as the **Path Conditions** as it tightens the bounds for paths chosen for regression; this is achieved by comparing the regression values of all viable paths and paths for which "far" in or out-of-the-money are discarded. The regression uses the paths between strike times for prediction. The issue for American options is that there is no "path" between the exercise times since every time is an expiration time. Thus, I decided that the data was too biased toward the most recent iterations to predict the next price accurately. 

A potential workaround is to use the $k$ previous iterations where $k$ is some user-selected number of iterations. I did not implement this as I attempted to keep this method as mathematically sound as possible.

Another issue is that this stratification requires access to **all** paths used for the calculation. [Loop blocking](https://www.intel.com/content/www/us/en/developer/articles/technical/efficient-use-of-tiling.html) was used to optimize runtime for the path generation so paths were generated in batches of around $200$, and accessing all paths was not possible.

### Loop Optimizations

[Loop blocking](https://www.intel.com/content/www/us/en/developer/articles/technical/efficient-use-of-tiling.html) (or tiling) is implemented when iterating over the total number of desired simulations by dividing the number of total simulations by $200$ to optimize cache allocation and hence speed. This is faster than a standard **for** loop since the blocks are small enough to be stored in CPU memory instead of DRAM, decreasing the "distance" the information has to travel to be accessed. The optimal block size is machine-dependent so I encourage any future users to test with their machine to find the most optimal block size. 

[Loop interchanges](https://en.wikipedia.org/wiki/Loop_interchange) were tested in the [*generate-data.cpp*](https://github.com/colalb1/Fast-Least-Squares-Monte-Carlo/blob/main/generate-data.cpp) file, but improvements were negligible. I assume this is because the order of the vectors being iterated over is in the single digits so changing them around would not have a great effect. The concept of "loop interchange" boils down to making the outermost **for** loop the loop with the greatest number of iterations/data and assigning the inner loop to that with fewer iterations/data since the fewer iterations will allow the data to stay in the CPU cache instead of moving to DRAM. This achieves a similar outcome to *loop blocking*. An example of this is the multiplication of a tall-and-skinny matrix $A$ by a square matrix $B$ as one desires the outer loop to iterate over $A$'s rows and the inner loop to iterate over $A$'s columns to keep the data in the CPU cache.

I also learned the term [loop fusion](https://en.wikipedia.org/wiki/Loop_fission_and_fusion) for which independent computations are put into the same **for** loop to save computation time. This is something I already do, but it's good to know the code is more efficient this way.

***************PARELLELIZATION*************** talk about thread pooling, mutex locking, computation time decrease


## Results

ADD METRICS AND PLOTS. In short, the optimized version that I wrote is nearly always more accurate than the original version. Regarding which basis is most accurate, it depends on the input parameters (ADD INFO TO BACK UP). In practice, one would likely choose the Power basis as it is easy to implement and provides marginally different results than the Hermitian and Laguerre bases. If you have an aching desire to know exactly which scenarios the Hermitian and Laguerre bases are more accurate, I encourage you to copy this code and proceed with rigorous testing.

## Conclusion

Add metrics on performance and discuss improvements in general.

## My Review of C++
Fast. Really slow to debug.
