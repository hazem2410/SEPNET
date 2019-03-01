<img src="images/logo.png" width = "500">

# Simulation and Estimation of Production Network

SEPNET is a tool to estimate the formation of links in a production network. The program is a simulation algorithm of the well known exponential random graph model [ERGM](https://www.sciencedirect.com/science/article/pii/S0378873306000372).
The code is based on the fixed-density Markov Chain Monte Carlo sampling and the stochastic approximation of the MLE. Details were discussed in [Snijders et al. (2006)](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1467-9531.2006.00176.x)
and [Morris et al. (2008)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2481518/). 

SEPNET implements ERGM for production network:

* SEPNET considers the gravity model assumption: Firms are connected based on their geographic distance and their respective performance (Sales, Profit, Productivity, etc.).
* SEPNET adds attributes related to the economic activity of firms (Industrial sector, Major bank, etc.)
* SEPNET is a flexible tool to add new economic statistics (Ownership relations, Bankruptcy, Financial Risk, etc.).
* SEPNET includes the most important network statistics (K- Triangles, K- stars, K- Two paths, etc.).

# 1. Usage details of SEPNET

* Download files: .cpp and .DAT for data.

* The program uses the C++ Eigen library.

* Eigen should be installed in your framework, e.g. directory C:/Eigen.

* To compile from your shell or cmd command: g++ -I C:/Eigen ERGMGithub.cpp -o ERGM.

* Then, simply execute ERGM.exe file.

# 2. Initialisation of parameters

Parameters in SEPNET considers mainly the chain size of the MCMC. Following, we describe each parameter with its suggested value.

* number_variables: It counts the number of attributes; 16 in our case.

* theta: The vector of estimated parameters, with size equal to number_variables.

* IFD_size: the size the fixed-density MCMC sampling. It depends on the network size, see [Lusher et al. (2013)](https://www.cambridge.org/core/books/exponential-random-graph-models-for-social-networks/9296EE2B53CDEF9FE9E2E981E2FDB8A8).

* a_r is the gain factor that controls the largeness of the updating steps of parameter theta. [Lusher et al. (2013)](https://www.cambridge.org/core/books/exponential-random-graph-models-for-social-networks/9296EE2B53CDEF9FE9E2E981E2FDB8A8)
suggests a value of 0.1.

* initial_simulation: the size of iteration to calculate the initial theta. This parameter is suggested to be equal to 3 + 7*number_variables.

* burnin: It defines the size of the burn-in stage.

* number_phases: The maximum number of iterations before convergence.

# 3. Publications based on SEPNET:

Krichene, H., Chakraborty, A., Fujiwara, Y., Inoue, H. and Terai, M. Tie-formation process within the communities of the Japanese production network: application of an exponential random graph model. Applied Network Science 4(5), (2019). https://doi.org/10.1007/s41109-019-0112-9 

Krichene, H., Arata, Y., Chakraborty, A., Fujiwara, Y. and Inoue, H. How Firms Choose their Partners in the Japanese Supplier-Customer Network? An application of the exponential random graph model (2018). RIETI DP 18-E-011. https://www.rieti.go.jp/jp/publications/dp/18e011.pdf