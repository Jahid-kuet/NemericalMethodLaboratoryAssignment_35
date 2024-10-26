# NemericalMethodLaboratoryAssignment_35
We have tried to create a console application implementing various neumeric methods . We have included Solution of Linear Equations, Solution of Non-linear Equations, Solution of Differential Equations and Matrix Inversion . We gave the user to choose various neumeric methods to solve linear equations such as  Jacobi iterative method,Gauss-Seidel iterative method, Gauss elimination ,Gausss-Jordan elimination,LU factorization. To solve non linear equations we gave the options Bi-section method, False position method,Secant method,
Newton-Raphson method. We have used Runge-Kutta method to solve differential Equations. We also gave the option to inverse a matrix .

In Jacobi iterative method,Gauss-Seidel iterative method  -> we have to carefully choose the augmented matrix so that it can be a "DIAGONALLY DOMINANT " matrix. 
In Gauss Elimination -> We performd row operations to transform the matrix into an upper triangular form, then back-substitutes to find solutions.
In Gauss-Jordan Elimination -> We did  further by transforming the matrix into reduced row-echelon form, making each variable's coefficient 1 and eliminating other entries in each column. Both functions print the solutions to the system.
In the lu method -> firstly we have created the lower and upper triangular matrix . Then we did the forward substitution to get the Y matrix and after doing backward substitution we got the solution matrix X.

<hr>

<h2>1) Matrix Inversion Using Gauss-Elimination and Gauss-Jordan Methods</h2>
<p><strong>Input format:</strong></p>
<pre>
Enter the number of unknowns (variables): n
Enter the augmented matrix (A | I form): A | I
</pre>

<p><strong>Output format:</strong></p>
<pre>
I | A<sup>-1</sup> 
where A<sup>-1</sup> is the inverse matrix of A.
</pre>

<h3>Test Case 1:</h3>
<p><strong>Input:</strong></p>
<pre>
3
1   2   3   1   0   0
4   6   7   0   1   0
-2  1   5   0   0   1
</pre>

<strong>Output:</strong>
<pre>
7.67    -2.33   -1.33
-11.33  3.67    1.67
5.33    -1.67   -0.67
</pre>

<h3>Test Case 2:</h3>
<strong>Input:</strong>
<pre>
4
0   0   -1  2   1   0   0   0
0   1   0   0   0   1   0   0
9   0   0   0   0   0   1   0
0   0   0   1   0   0   0   1
</pre>

<p><strong>Output:</strong></p>
<pre>
0.00    0.00    0.11    0.00
0.00    1.00    0.00    0.00
-1.00   0.00    0.00    2.00
0.00    0.00    0.00    1.00
</pre>

<hr>

<p><strong>2) Solution of Differential Equations:</strong></p>
<p><strong>Runge-Kutta Method (RK4) </strong></p>
Consider the Differential Equation: dy/dx = x - y

Input:
Enter h(step size) and x(final value of x):
h = 1  
x = 10

Output:
Tabular representation of X, Y (RK4), Y (Exact), and Error
Bisection method :  In this method, the interval distance between the initial values is treated as a line segment. It then successively divides the interval in half and replaces one endpoint with the midpoint so that the root is bracketed. This method is based on The Intermediate Value Theorem.

False Position method: Step 1: we choosed two initial points a and b such that function at those points have opposite sign  f(a)⋅f(b)<0.

Step 2: we calculated the point c where the linear approximation intersects the x-axis using the formula.

Step 3: we determined f(c).

i) If f(c) ⋅ f(a) < 0, then the root lies between a and c. Set b = c.
ii) If f(c) ⋅ f(b) < 0, then the root lies between b and c. Set a = c.
Step 4: we repeatd the steps until ∣f(c)∣ is less than a previously defined tolerance level .

Secant method:
1. We started with two initial guesses, x0 and x1 to find the root.
2. We calculated a new approximation by the secant formula.
3. We updated guesses until the difference between successive approximations is within a given  tolerance.
4. repeatd until the root is reached.

Newton-Raphson method: 
1. We started from an initial guess.
2. We updated it using the formula xnew = x - f(x) / f'(x)
3.  set value xcurrent=xnew.
4. Continue applying the formula until the difference between successive approximations is within a chosen tolerance level.







