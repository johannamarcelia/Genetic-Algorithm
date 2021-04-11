# Genetic-Algorithm
This contains the functions required to run the genetic algorithm for changepoint detection. These functions were used in the analysis for "A Data Adaptive Model for Retail Sales of Electricity."

The function takes in time series data, and was specifically constructed to work with the EIA retail sales of electricity to consumer by end use sector. The first function computes the log liklihood values of the selected regression model while the second function calculates the corresponding MDL penalty. The third function combines these into a "fitness function", and the final function contains the genetic algorithm. 


