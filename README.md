# Genetic-Algorithm
This contains the functions required to run the genetic algorithm for changepoint detection. These functions were used in the analysis for "A Data Adaptive Model for Retail Sales of Electricity."

The function takes in time series data, and was specifically constructed to work with the EIA retail sales of electricity to consumer by end use sector. The first function computes the log liklihood values of the selected regression model while the second function calculates the corresponding MDL penalty. The third function combines these into a "fitness function", and the final function contains the genetic algorithm. 



Data fitting function: 
Fits model to the data for some given changepoint configuration, and provides a likelihood value that will be used to assess the fit.
Inputs: x = time series formatted dataset
        crm = a "chromosome" of the number of changepoints and their locations in the data set
Outputs: dfv = the logliklihood value, and the residual values that pertain to that model


Penalty function: 
Computes penalty term based on the number of parameters requried to fit model for a given chromosome configuration. Based on MDL theory
Inputs: x = time series formatted dataset
        crm = a "chromosome" of the number of changepoints and their locations in the data set
Outputs: pnt = pentalty value


Fitness function: 
This function will be optimized when calculated with the best possible configuration of changepoints. Combines the data fitting and penalty functions to arrive at a single "fitness score" value
Inputs: x = time series formatted dataset
        crm = a "chromosome" of the number of changepoints and their locations in the data set
Outputs: the penalty value for a given chromosome configuration


Changepoint function: 
The genetic algorithm is used to stocastically test many different changepoint configurations and iteratively arrive at the configuration that best optimizes the fitness function 
Inputs: x = time series formatted data
Outputs: params = model intercept terms, model slope terms, optimal changepoint configuration


Known Issues: for R to execute the arima function in line 41, scaling of the data is required. Additionally, the method used greately impacts the speed with which the algorithm runs. CSS is the fastest, but if a warning or error occures, the data may need to be scaled further, or a different method (CSS-ML) can be used. 


