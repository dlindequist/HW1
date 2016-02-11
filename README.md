# HW1

- the main code is in Hw1.m 
- to run the code you need root_MP.m which contains the fsolve routine

Main idea of the coding: I use Tauchen's method to construct the grid for z (50 grid points) and the Markov matrix P. I also programmed 
Rouwenhorst's method, but I am not using it here. Then I use an fsolve routine to find theta's corresponding to the z's solving the MP model.
In a next step, I simulate 50 observations of the productivity shock from the Markov chain. These shocks are then linked to their 
corresponding thetas. Finally, I calculate job finding and vacancy filling probabilities as well as implied unemployment rates.
To compare the baseline specification with the Hagedorn Manovskii specification, change the parameter values in line 58 and 60 of the code.
Plots are attached. 
