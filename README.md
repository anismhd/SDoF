# Single Degree of Freedom System Solver 

\begin{equation}
	M \ddot{x} + C \dot{x} + K x = F(t)
\end{equation}

## Time Integration or Frequency Domain Analysis 
This is function written in python to solve linear single degree of freedom system
The theoretical background of the same can be found in following book

	1 - An Introduction To The Finite Element Method,*McGraw-Hill series in mechanical engineering*,J. N. Reddy, Edition 3, 2006, ISBN	007053084X, 9780070530843, 766 pages

| Inputs  | Description             |
| ------- | ----------------------- |
| *N*     | define type of analysis |
|         | 0 - the constant-average acceleration method (stable) |
|         | 1 - the linear acceleration method (conditionally stable) |
|         | 2 - the central difference method (conditionally stable) - This shouldn't be used |
|         | 3 - the Galerkin method (stable) |
|         | 4 - the backward difference method (stable) |
|*M*   | single degree of freedom mass |
|*C*| single degree of freedom damping co-efficient |
|*K*| stiffness of single degree of freedom |
|*t*| time index for forcing time history |
|*F*| Forcing time history |
|$x_0$ | Initial displacement |
|$v_0$ | Initial velocity |

Inputs		:: N, Mass, Damping, Stiffness, Time, Acceleration, Initial displacement, Initial velocity
Where N will be used to define type of analysis
	0 - the constant-average accelaration method (stable)
	1 - the linear accelaration method (conditionally stable)
	2 - the central difference method (conditionally stable) - This should'nt be used
	3 - the Galerkin method (stable)
	4 - the backward difference method (stable)

Output		:: displacement, velocity, accelaration