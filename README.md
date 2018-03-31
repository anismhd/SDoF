# Single Degree of Freedom System Solver 

Anis Mohammed Vengasseri

anis.mhd@gmail.com

## Time Integration or Frequency Domain Analysis 
This is function written in python to solve linear single degree of freedom system
The theoretical backhround of the same can be found in following book

Title	An Introduction To The Finite Element Method
McGraw-Hill series in mechanical engineering
Author	J. N. Reddy
Edition	3
Publisher	McGraw-Hill, 2006
ISBN	007053084X, 9780070530843
Length	766 pages

	|Inputs | Description|
	|-------|------------|
	|*N*| define type of analysis|

	Markdown | Less | Pretty
	--- | --- | ---
	*Still* | `renders` | **nicely**
	1 | 2 | 3

Inputs		:: N, Mass, Damping, Stiffness, Time, Accelaration, Initial displacement, Initial velocity
Where N will be used to define type of analysis
	0 - the constant-average accelaration method (stable)
	1 - the linear accelaration method (conditionally stable)
	2 - the central difference method (conditionally stable) - This should'nt be used
	3 - the Galerkin method (stable)
	4 - the backward difference method (stable)

Output		:: displacement, velocity, accelaration