# Fixed Point Updater in MATLAB

**Author:** Alex Clymo  

This repository provides a flexible and modular toolkit for solving fixed point problems of the form `x = f(x)` in MATLAB, where `x` is a column vector of length `N`. The methods are all iterative, using the current guess `x_k` and evaluation `f(x_k)` to build the new guess `x_{k+1}`, allowing the code to be implemented in a simple loop. The code implements adaptive dampening, and certain methods automatically store a history of past guesses and evaluations in order to accelerate convergence by, for example, approximating the Jacobian. 

It can also be applied to solving nonlinear equations of the form `g(x) = 0` by simply adding `x` or `-x` to both sides and therefore defining `f(x) = x + g(x)` or `f(x) = x - g(x)`. For Jacobian based methods, either definition is fine, while for the basic dampened fixed point update which version you choose matters.

> 🚧 **Warning!** This code is very much in early development. I put it online at this early stage to encourage myself to start using Github. Please use with caution, and comments are always welcome. 

## 🔧 Overview

Solving fixed point problems is a core task in many quantitative macroeconomics applications. This toolkit is especially useful in the following scenarios:

- Embedding fixed point iterations within a broader calibration or simulation framework. E.g.
    - Wrapping a calibration loop around the solution of the steady state of your model
    - Solving for the equilibrium price sequence following an MIT shock
- When putting your model into a Matlab function to use `fsolve` is inconvenient.
- Testing and comparing update methods like damped iteration or Anderson acceleration.

The code is not particularly sophisticated, but the idea is for it to be very practically useful. If you have an updating scheme where you currently update some vector `x` in a loop with dampening, this can be slow. This function is meant to allow you to replace that dampened update with something more sophisticated with little to no hassle. 

## 📁 File Descriptions

| File               | Description |
|--------------------|-------------|
| `fpUpdate.m`       | Core function for updating the guess `x` using various fixed point methods. Manages internal histories needed for acceleration and Jacobian methods. |
| `fpSetup.m`        | Helper function to initialize the `par` structure with default parameters for the chosen method. |
| `example1_basics.m`| Script demonstrating how to use the fixed point solver with example functions, comparing performance of methods and validating against MATLAB’s `fsolve`. |

## 🚀 Usage

The core use is to update a guess `x` to a new guess `x_new` using the function call:
```matlab
[x_new, par] = fpUpdate(x, fx, par);
```
where `fx` is the value of `f(x)` evaluated at `x` (i.e. `fx = f(x)`). The structure `par` defines the update method, contains the method's options, and automatically stores and updates any past function evaluation values or data needed to perform the method. Note that `par` is updated with new data as part of the function's output.

To use this code, the `par` structure must be created towards the top of your code. This is done using the `fpSetup` function. For example, to set up the solver to use Anderson Acceleration, we automatically set up `par` with the default parameters using
```matlab
par = fpSetup('anderson');
```
The parameters can then be manually edited if desired by editing the `par` structure. 

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
As stated above, this code is particularly useful when you do not want to put the `f(x)` function into a Matlab function. If you are happy to save your problem into a function `f` so that `fx = f(x);` you are probably better off doing so and sending the function to `fsolve`. 
But for many practical applications you might want to keep your model code in the same script as the update step (at least during code development). In this case, `fpUpdate` might be useful, and can be faster than running `fsolve` out of the box on a function with thousands of variables (such as a price sequence).

## 📈 Methods Currently Supported

All methods optionally implement adaptive dampening, where the dampening parameter `zeta` is lowered (raised) if the error is rising (falling).

- **Fixed Point with dampening**:
    - Simple update `x_new = zeta .* fx + (1 - zeta) .* x` where `zeta` is a dampening parameter which can be either a scalar or vector of length `N`.
    - Works for contractions, possibly fails if not. 
    - In price sequence update loops, often leads to oscillations and overshoots.
- **Anderson Acceleration**:
    - Uses history of past guesses and residuals to improve convergence. Solves a least squares problem each iteration to choose weights and updates `x_new` as a weighted sum of past function evaluations. See [here](https://en.wikipedia.org/wiki/Anderson_acceleration) for details.
    - Smaller memory requirement than Jacobian based methods, while still improving speed. Idea is that the partial history approximates the role of the Jacobian. For example, in price sequence update loops, this extra information avoids oscillations and overshoots, and allows you to use less dampening than the simple fixed point method and so converge in fewer iterations.
    - Code automatically stores history of last `Ma` guesses and and function evaluations in `par`. 
    - Seems typical to set `Ma` to around 5 or 10. 
- **Broyden's Method**:
    - Jacobian based quasi-Newton method: builds an approximation to the inverse Jacobian using the history of past guesses and function evaluations. 
    - Higher memory requirement than fixed point or Anderson method when `N` is large. Might be infeasible for, e.g., solving long price sequences.
    - Code automatically stores and updates inverse Jacobian in `par`. Uses Broyden's first (i.e. good) method, and works directly with inverse Jacobian. See [here](https://en.wikipedia.org/wiki/Broyden%27s_method) for details.
    - Code includes protection against small denominator in Jacobian update, and automatically resets Jacobian if it becomes ill-conditioned.
- **Diagonal Jacobian**
    - Method assumed Jacobian is diagonal and simply estimates it using a finite difference approximation, comparing current evaluation to past evaluation. 

## ⚙️ Requirements

- MATLAB (tested with R2024b)
- Optimization Toolbox (optional, for `lsqlin` used in Anderson acceleration)
