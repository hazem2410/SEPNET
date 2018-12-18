# Simulation and Estimation of Production Network

SEPNET is a tool to estimate the formation of links in a production network. The program is a simulation algorithm of the well known exponential random graph model [ERGM](https://www.sciencedirect.com/science/article/pii/S0378873306000372).
The code is based on the fixed-density Markov Chain Monte Carlo sampling and the stochastic approximation of the MLE. Details were discussed in [Snijders et al. (2006)](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1467-9531.2006.00176.x)
and [Morris et al. (2008)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2481518/). 

SEPNET implements ERGM adapted for production network:

* SEPNET considers the gravity model assumption: Firms are connected based on their geographic distance and their respective performance (Sales, Profit, Productivity, etc.).
* SEPNET adds attributes related to the economic activity of firms (Industrial sector, Major bank, etc.)
* SEPNET is a flexible tool to add new economic statistics (Ownership relations, Bankruptcy, Financial Risk, etc.).
* SEPNET includes the most important network statistics (K- Triangles, K- stars, K- Two paths, etc.).