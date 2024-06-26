# MatrixStabilization
Given an unstable matrix $A$, i.e. with some eigenvalues in the right complex half-plane, it computes an approximation of the stable matrix $B=A+\Delta$ such that $\||\Delta\||_F=1$ is minimized. 

The codes implemented here refers to the paper

N. Guglielmi, S.Sicilia, "Stabilization of a matrix via a low-rank-adaptive ODE", https://arxiv.org/html/2402.14657v1

The main function that returns the distance and the perturbation is 'Stabilize' and it can solve the both the unstructured and strucured problem, with adaptive-rank or fixed-rank integrator.

The folder 'utils' contains all the subfunctions needed by the main function 'Stabilize'.

The folder 'tests' contains the numerical experiments of the paper. 

Before running any code, please run Install.m to have all paths added. 



