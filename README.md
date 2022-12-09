# Regime-switching-in-commodity-prices
This repo contains the codes of [Wars, cartels and COVID-19: regime switching in commodity prices](https://doi.org/10.1080/13504851.2022.2133892). 

We follow [Hamilton (1989)](https://doi.org/10.2307/1912559 ) and propose a tractable approach to modelling changes in regime as the outcome of a discrete-state Markov process. The mean and the variance of the real price of oil and copper are parameterized in terms of an unobserved state variable that follows a 3-state Markov process with unknown transition probabilities. We want to determine the probability that the variance of either, the price of oil or copper is, at any point in time, in a given state.

 **Description of files** 
- `Gen_TSMS.jl`: generate the structure
- `TSMS_functions.jl`: auxiliary functions
- `POIL_PCU_2021.csv`: the data
- `3States_Markov_Switching-OIL-CUP.ipynb`: Jupyter notebook that describes and runs all the code
