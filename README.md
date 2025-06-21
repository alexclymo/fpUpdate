# Fixed Point Updater in MATLAB

**Author:** Alex Clymo  
**Date:** 18 June 2025

This repository provides a flexible and modular toolkit for solving fixed point problems of the form $x = f(x)$ in MATLAB, where $x$ is a column vector of length $N$. It can also be applied to solving nonlinear equations of the form $g(x) = 0$ by simply defining $f(x) = g(x) + x$. It includes a generic solver wrapper `fpUpdate`, setup utility `fpSetup`, and an example script demonstrating their usage.

> 🚧 **Warning!** This code is very much in early development. I put it online at this early stage to encourage myself to start using Github. Please use with caution, and comments are always welcome. 

## 🔧 Overview

Solving fixed point problems is a core task in many quantitative macroeconomics applications. This toolkit is especially useful in the following scenarios:

- Embedding fixed point iterations within a broader calibration or simulation framework. E.g.
    - Wrapping a calibration loop around the solution of the steady state of your model
    - Solving for the equilibrium price sequence following an MIT shock
- When putting your model into a Matlab function so you can use `fsolve` is inconvenient.
- Testing and comparing update methods like damped iteration or Anderson acceleration.

The code is not particularly sophisticated, but the idea is for it to be very practically useful. If you have an updating scheme where you currently update some vector $x$ in a loop with dampening, this can be slow. This function is meant to allow you to replace that dampened update with something more sophisticated with little to no hassle. 

## 📁 File Descriptions

| File               | Description |
|--------------------|-------------|
| `fpUpdate.m`       | Core function for updating the guess `x` using various fixed point methods (`fixedPoint`, `anderson`). Manages internal histories needed for acceleration methods. |
| `fpSetup.m`        | Helper function to initialize the `par` structure with default parameters for the chosen method. |
| `example1_basics.m`| Script demonstrating how to use the fixed point solver with example functions, comparing performance of methods and validating against MATLAB’s `fsolve`. |

## 🚀 Usage

The core use is to update a guess `x` to a new guess `x_new` using the function call:
```matlab
[x_new, par] = fpUpdate(x, fx, par);
```
where `fx` is the value of $f(x)$ evaluated at `x` (i.e. `fx = f(x)`). The structure `par` defines the update method, contains the method's options, and automatically stores and updates any past function evaluation values or data needed to perform the method. Note that `par` is updated with new data as part of the function's output.

To use this code, the `par` structure must be created towards the top of your code. This can be done either manually or using the `fpSetup` function. For example, to set up the solver to use Anderson Acceleration, we could manually type:
```matlab
par.method = 'anderson';
par.Ma = 5; %number of last guesses to use in Anderson scheme
par.zeta0 = 0.01; %dampening during pre-Anderson phase
par.zeta1 = 1; %dampening during Anderson phase
```
Or automatically set up `par` with the default parameters using
```matlab
par = fpSetup('anderson');
```
The parameters can then be manually edited if desired.

To use the solver, it is then placed into a standard `while` loop where we update `x` until convergence:
```matlab

% initial guess
x0 = ...
x = x0;

diff = 1;
tol = 1e-5;
while diff > tol

    % evaluate function at current guess
    fx = ...

    % Perform one x update to get x_new
    [x_new,par] = fpUpdate(x,fx,par);

    % compute percentage error (add small constant to denom to handle x=0)
    diff = max( abs( (y-x)./(1e-3+abs(x)) ) )

    % update x by setting x = x_new (only if not converged yet)
    if diff > tol
        x = x_new;
    end

end
```
As stated above, this code is particularly useful when you do not want to put the $f(x)$ function into a Matlab function. If you are happy to save your problem into a function `f` so that `fx = f(x);` you are probably better off doing so and sending the function to `fsolve`. 
But for many practical applications you might want to keep your model code in the same script as the update step (at least during code development). In this case, `fpUpdate` might be useful, and can be faster than running `fsolve` out of the box on a function with thousands of variables (such as a price sequence).

## 📈 Methods Currently Supported

- **Fixed Point with dampening**: $x_{k+1} = \zeta \odot f(x_k) + (1-\zeta) \odot x_k$
    - Simple update `x_new = zeta .* fx + (1 - zeta) .* x` where `zeta` is a dampening parameter which can be either a scalar or vector of length `N`.
    - Works for contractions, possibly fails if not. 
    - In price sequence update loops, often leads to oscillations and overshoots.
- **Anderson Acceleration**: $`x_{k+1} = \sum_{i=0}^{m} (\alpha_k)_i f_{k - m + i}`$
    - Uses history of past guesses and residuals to improve convergence. Solves a least squares problem each iteration to choose parameters $\alpha_k$ and updates `x_new` as a weighted sum of past function evaluations. See [here](https://en.wikipedia.org/wiki/Anderson_acceleration) for details.
    - Smaller memory requirement than Jacobian based methods, while still improving speed. Idea is that the partial history approximates the role of the Jacobian. For example, in price sequence update loops, this extra information avoids oscillations and overshoots, and allows you to use less dampening than the simple fixed point method and so converge in fewer iterations.
    - Code automatically stores history of last `Ma` guesses and and function evaluations in `par`. Code allows for standard dampening on top of the Anderson update. 
    - Seems typical to set `Ma` to around 5 or 10. 
    - Convergence is not guaranteed, and after getting close to the solution quite quickly, the method might get stuck at, e.g., an error of around 1e-3. 
- **Broyden's Method (in progress!)**:
    - Jacobian based method quasi-Newton method: builds an approximation to the inverse Jacobian using the history of past guesses and function evaluations. 
    - Higher memory requirement than fixed point or Anderson method when $N$ is large. Might be infeasible for, e.g., solving long price sequences.
    - Code automatically stores and updates inverse Jacobian in `par`.
- **Limited-memory BFGS (L-BFGS) Method (in progress!)**
    - A limited memory version of Broyden-type methods, which only stores a partial history of function evaluations.
    - Memory requireement therefore same as Anderson Acceleration, and might be useful in situations where $N$ is large and Broyden is infeasible.
    - Code automatically stores history of last guesses and and function evaluations in `par`.

## ⚙️ Requirements

- MATLAB (tested with R2024b)
- Optimization Toolbox (optional, for `lsqlin` used in Anderson acceleration)



## 💡 Regularisation in Anderson Acceleration

When the residual history matrix $`R`$ in Anderson Acceleration becomes ill-conditioned, the least squares step can become numerically unstable. This happens as we approach the solution, and all the columns of $R$ become very similar, amplifying the role of any noise in the calculation of the $\alpha$ weights. To improve stability, we apply **Tikhonov regularisation** (ridge regression) to the Anderson step.

### Regularised Anderson Problem

We solve the following constrained optimisation problem:

$$\min_\alpha ||R \alpha||^2 + \lambda^2 ||\alpha||^2 
\text{ subject to } \sum_i \alpha_i = 1$$

The addition of the $\lambda^2 ||\alpha||^2$ term stabilises the solution when columns of $`R`$ are nearly linearly dependent. When $\lambda=0$ this returns to the standard problem, and when $\lambda\rightarrow\infty$ the solution approaches equal weights ($`\alpha = (1,1,...,1)/\text{length}(\alpha)`$).

At iteration $k$, with the solved $\alpha_k$ vector in hand, the Anderson method updates to the next guess $x_{k+1}$ using the formula $`x_{k+1} = \sum_{i=0}^{m} (\alpha_k)_i f_{k - m + i}`$

### Condition Number of the Regularised Matrix

Let $`\sigma_1 \geq \dots \geq \sigma_n`$ be the singular values of $`R`$. The condition number of the matrix $`R'R + \lambda^2 I`$ is:

$$\kappa = (\sigma_1^2 + \lambda^2) / (\sigma_n^2 + \lambda^2)$$

- The condition number of a matrix equals 1 if it is perfectly conditioned, and approaches infinity the columns become perfectly collinear.
- The condition number of $R'R$ is equal to the condition number of $R$ squared.
- As $`\lambda \to 0`$, this approaches the (squared) condition number of $`R`$. As $`\lambda \to \infty`$, $`\kappa \to 1`$.

### Choosing $`\lambda`$ to Target a Desired Condition Number

To choose $`\lambda`$ such that the condition number of the regularised matrix equals a target value $`\kappa_{\text{target}}`$, use:

$$\lambda^2 = (\sigma_1^2 - \kappa_{\text{target}} \cdot \sigma_n^2) / (\kappa_{\text{target}} - 1)$$

We impose a maximum value of the condition number of $R$ we will tolerate, `par.maxCondR`. If $R$ is worse conditioned than this, then $\lambda$ is chosen according to the above formula. If not, then $\lambda=0$ and no regularisation is applied. As `par.maxCondR` approaches infinity, regularisation is turned off. As `par.maxCondR` approaches 1, $\alpha$ is forced to be a vector of equal weights.
