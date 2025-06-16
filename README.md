# HSA-DC
The hybrid simulated annealing approach for global minimization of DC functions

HSA-DC is a hybrid optimization method (Fortran 77) coded by Kaisa Joki and Adil M. Bagirov to solve global minimization problems where the objective function is expressed as a difference of two convex (DC) functions, subject to box constraints. The method integrates a local DC optimization technique with a simulated annealing (SA) mechanism. Initially, it identifies a critical point using a local search. Then, it employs approximated epsilon-subdifferentials of the DC components —- represented by polytopes —- to explore candidate points with potentially lower objective function values. Since these approximations may not always lead to global descent directions, the method incorporates the SA acceptance criterion to probabilistically accept or reject new points, enhancing the chances of escaping local minima while maintaining computational efficiency.

The HSA-DC method utilizes the aggregate subgradient method (AGGSUB) from [2] as a local search method. Additionally, the augmented subgradient method (AUGSUB) [3] is employed during the last iteration round to further improve the solution. To solve quadratic programming problems, the algorithm from [4] is used.

The software is free for academic teaching and research purposes but we ask you to refer the reference given below if you use it. To use the software modify HSA-DC_method.for as needed. If you have any questions conserning the software, please contact directly the author Kaisa Joki (email: kjjoki@utu.fi).

# Code include:                                                                     
         
   HSA-DC_method.for    - Main program for the HSA-DC method           
                                                                                              
   Makefile             - Makefile       

   startpoints5.txt     - Random starting points from the boxes [-5,5]
   
   startpoints10.txt    - Random starting points from the boxes [-10,10]
   
   startpoints100.txt   - Random starting points from the boxes [-100,100]
   
   
# References:                                                                        
                                                                                              
[1] Bagirov, A.M., Joki, K., Mäkelä, M.M, Taheri, S. "A hybrid simulated annealing approach for global minimization of DC functions". (manuscript)

[2] Bagirov, A.M., Taheri, S., Joki, K., Karmitsa, N., Mäkelä, M.M.: "Aggregate subgradient method for nonsmooth DC optimization". Optim. Lett. 15(1), pp. 83–96 (2021).                                                   

[3] Bagirov, A.M., Hoseini Monjezi, N., Taheri, S.: "An augmented subgradient method for minimizing nonsmooth DC functions". Comput. Optim. Appl. 80(1), pp. 411–438 (2021).                                               

[4] Wolfe, P.H.: "Finding the nearest point in a polytope". Math. Program. 11(2), pp. 128–149 (1976). 
