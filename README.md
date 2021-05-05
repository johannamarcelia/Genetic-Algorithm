# Genetic-Algorithm
This contains the functions required to run the genetic algorithm for changepoint detection. These functions were used in the analysis for "A Data Adaptive Model for Retail Sales of Electricity."

The function takes in time series data, and was specifically constructed to work with the EIA retail sales of electricity to consumer by end use sector. The first function computes the log liklihood values of the selected regression model while the second function calculates the corresponding MDL penalty. The third function combines these into a "fitness function", and the final function contains the genetic algorithm. 

Inputs: Time series formatted data
Outputs:

Known Issues: for R to execute the arima function in line XX, scaling of the data is required. Additionally, the method used greately impacts the speed with which the algorithm runs. CSS is the fastest, but if a warning or error occures, the data may need to be scaled further, or a different method (CSS-ML) can be used. 


