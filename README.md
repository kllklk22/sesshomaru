# sesshomaru

high-dimensional pde solver for option pricing. wrote this at 4am after watching my professor's python code take ten minutes to price a simple european derivative.

this uses an implicit finite difference method to solve the black-scholes equation: 
$\frac{\partial V}{\partial t} + \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 V}{\partial S^2} + rS\frac{\partial V}{\partial S} - rV = 0$

it's vectorized and uses julia's native linear algebra backend. it prices grids instantly. 
