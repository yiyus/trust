# Trust-Region Newton Minimisation

This is a trust-region Newton minimisation framework for Dyalog APL. When the Hessian matrix is not available, an estimation given by the Broyden–Fletcher–Goldfarb–Shanno (BFGS) or Levenberg-Marquardt (LM) algorithms is used, depending on whether the problem consists of minimising a scalar function or a vector of residuals (with an associated cost function).

The package is intended to be robust and offer maximum flexibility, while at the same time providing maximum convenience to the user.

## Introduction

The goal of the Newton method is to find the set of *solution parameters* $p$ which minimises an objective (or cost) function $f(p)$. Each iteration solves the damped linear system

$$(H + \lambda D)\Delta p = -g$$

for a step $\Delta p$, where $g$ is the gradient of $f$ at $p$, $H$ is the (positive-semidefinite) Hessian, $\lambda$ is a damping factor and $D$ is a diagonal scaling matrix with $D_{kk} = \max(\epsilon(\lambda), H_{kk})$ and $\epsilon(\lambda)$ an adaptive damping floor. The damping factor $\lambda$ interpolates between the aggressive Newton-like step (small $\lambda$) and a small gradient-descent step (large $\lambda$), defining an implicit trust region.

For the new guess $p + \Delta p$, the predicted error reduction is

$$\Delta s_p = -\tfrac{1}{2}(\Delta p^T g - \lambda \Delta p^T D \Delta p)$$

If the ratio of the actual error reduction obtained at the next iteration with respect to this prediction is above a configurable gain threshold, the guess is accepted and $\lambda$ is decreased; otherwise the guess is rejected and $\lambda$ is increased.

The algorithm finishes once one of the following conditions is met, returning the last accepted guess:

- Maximum number of iterations reached
- Objective function below specified tolerance
- Relative change (in parameters or objective function) below specified tolerance
- Stagnation at maximum damping factor

### Hessian and gradient estimation

When the true Hessian is not available, it can be estimated applying the BFGS or LM methods.

The BFGS method provides an approximation of the Hessian when the objective is to minimise a scalar function. For the accepted solution $k$, the Hessian is estimated as the matrix $H_k$ obtained such that:

$$H_{k+1}=H_k + \frac{\Delta g \Delta g^T}{\Delta g^T \Delta p} - \frac{H_k \Delta p \Delta p^T H_k^T}{\Delta p^T H_k \Delta p}$$

where $\Delta p$ and $\Delta g$ are the changes in parameters and gradient from step $k$ to step $k+1$, and starting with $H_0 = I$.

If the goal is to minimise residuals instead of a scalar function, the LM algorithm can be applied to minimise a cost function $\sum_i \rho(y_i(p), c_i)$, where $\rho$ is a loss function which depends on the residuals $y_i(p)$ and arbitrary weight factors $c_i$.

The Hessian and gradient of the cost function

$$H = J^T W J$$
$$g = J^T W y$$

depend on the Jacobian of the function to optimise, $J$, and a weight matrix $W$ which depends on the choice of loss function (the identity matrix for $L_2$).

#### Loss Functions

The standard loss function $L_2$ calculates the loss as the square of the residual, using the identity matrix for weights. Robust loss functions provide mechanisms to mitigate the effect of outliers in fitting data:

* **Huber** Equivalent to $L_2$ for small residuals and linear for larger residuals
* **Cauchy** Strongly down-weights large outliers
* **SoftL1** Smooth approximation of $L_1$ loss that behaves like $L_2$ for small residuals
* **Tukey** Redescending M-estimator that completely rejects extreme outliers
* **Welsch** Another redescending M-estimator, smoother than Tukey's in its rejection
* **Fair** Less sensitive to large errors than $L_2$, but not redescending
* **Arctan** Limits maximum loss of individual residuals

### Bounded parameters

Parameter bounds can be enforced via a trust-region interior-point method based on Coleman and Li's affine-scaling approach.

At each iteration, parameters are clamped to the interior of the feasible region with a small margin. A diagonal penalty matrix $\Delta H_{kk} = v_k^{-1/2}$ is added to the Hessian, where $v_k$ is the Coleman-Li distance (how far each parameter is from the bound it is heading towards):

$$v_k = \begin{cases} u_k - p_k & \text{if } g_k < 0 \\ p_k - l_k & \text{if } g_k > 0 \end{cases}$$

### Normalised damping factor

The damping factor $\lambda$ will oscillate between the specified limits $\lambda_{min}$ and $\lambda_{max}$, starting with an initial value $\lambda_0$. The *normalised damping factor* is calculated as:

$$\lambda_{norm} = \frac{(\lambda_{max}-\lambda_0)(\lambda-\lambda_{min})}{(\lambda_0-\lambda_{min})(\lambda_{max}-\lambda)}$$

If $\lambda_{norm}$ is 1, the initial damping factor is being used. If it is larger, it means that more damping than indicated was necessary; if it is lower, it means that damping could be decreased for faster convergence. The use of normalised factors between zero and infinity allows tracking and adjusting the effective size of the trust region independently of the particular quantities of the problem.

### Adaptive damping floor

To ensure stability when dealing with highly ill-conditioned or rank-deficient Hessian matrices, an adaptive floor $\epsilon(\lambda)$ is applied to the diagonal scaling matrix $D$.

The floor takes an arbitrarily small value $\epsilon_0$ during optimistic low-damping phases ($\lambda_{norm} \le 1$), but smoothly scales toward 1 as overall damping approaches $\lambda_{max}$. This ensures weak parameter components are strongly buffered when navigating difficult topological features.

## Usage

### `Eval` Operator

    R←{X}f Eval Y

`Eval` is a universal oracle operator that takes a monadic objective function `f` as its left operand and a set of parameters `Y`, and returns the corresponding residuals and derivative (either a gradient vector, or a Jacobian matrix). If `f Y` returns a nested vector of more than one element, it is returned as is. Otherwise, the result is treated as a vector of residuals and the Jacobian matrix (or gradient vector) is estimated using finite differences. The optional left argument `X` sets the relative perturbation step used for the numerical estimation (defaulting to `⎕CT*÷2`).

### `Min` Operator

    R←{X}f Min Y

`Min` is a monadic operator. It takes a left operand to return an ambivalent function which minimises an objective function given an initial set of parameters. Several configuration options are available, with sensible defaults previously defined.

The left operand `f` must be an evaluation function or a configuration namespace containing an `Eval` function. The return value of `f` determines the optimisation algorithm to use.

- If `f` returns a vector, they are interpreted as residuals and LM is used to minimise the cost function. `f` might also return the corresponding Jacobian, else it is estimated numerically.
- If `f` returns a scalar value, BFGS is used. The gradient is estimated numerically if not returned by `f`.
- If `f` returns a scalar value, Hessian matrix, and gradient, the Newton method is used.

`Y` must be a vector. The first element of `Y`, or `⊂Y` if `1=≡Y`, contains the initial guess for the solution parameters.
If the next element of `Y` is a scalar numeric value, it is interpreted as the initial normalised damping factor.
Additional elements of `Y` must be configuration namespaces. The final configuration parameters are obtained by overwriting the parameters in the namespace given as left operand with those given as right argument from right to left. Default values will be used for non-defined parameters, except for `Eval`, which must be provided either as left operand or as member of a configuration namespace.

Configuration namespaces may define the following options:

* `toli`: Maximum number of iterations (default `1E3`)
* `tolc`: Tolerance for convergence of objective function (default `⎕CT`)
* `tolr`: Tolerance for relative change, either in the solution parameters or the objective function (default `⎕CT`)
* `tolg`: Tolerance for the gain ratio to accept or reject a guess (default `1E¯2`)
* `dini`: Initial damping factor for `dnorm=1` (default `1E¯2`)
* `dinc`: Multiplier to increment damping factor after rejected guess (default `5`)
* `ddec`: Multiplier to decrement damping factor after accepted guess (default `÷dinc`)
* `dmax`: Maximum damping factor (default `÷⎕CT`)
* `dmin`: Minimum damping factor (default `÷dmax`)
* `pert`: Relative perturbation applied during finite-difference estimation of derivative (default `⎕CT*÷2`)
* `loss`: Loss function for vector residuals: `L2` `Huber` `Cauchy` `SoftL1` `Tukey` `Welsch` `Fair` `Arctan` or dyadic function (default `L2`)
* `scale`: Scale factor passed as left argument to loss function (default for 95% efficiency in robust loss functions)
* `verbose`: If `1`, print `iter cost rel dnorm p` each iteration (default `0`)

Configuration namespaces may also contain an optional `Callback` function executed prior to convergence checks.

The following options apply only to bounded problems (if no bounds are set, the solver runs unconstrained):

* `lower`: Lower bounds for the solution parameters.
* `upper`: Upper bounds for the solution parameters.
* `margin`: Margin from the bounds for the interior-point clamp.

The returned value `R` is a namespace including all the configuration options and the additional fields:

* `iter`: Total number of iterations
* `cost`: Final value of objective function
* `rel`: Final relative change metric
* `dnorm`: Final normalised damping factor
* `p0`: Initial guess of solution parameters
* `p`: Final optimised solution parameters

#### Notes

* With the exception of `toli` and `tolc`, configuration parameters should generally be modified only by expert users or in case of convergence problems

* The relative change metric `rel` is the minimum relative change between successive accepted solutions either in the cost or in the solution parameters

* In addition to being used for the definition of default values, `⎕CT` is also the baseline for adaptive floor damping

* The perturbation to estimate the Jacobian `pert` and the scaling factor for loss functions `scale` can be either scalar values, or vectors of the same length of respectively the parameters and the residuals

* Loss functions and their respective weights, as well as their default values for the scaling parameter are defined in the namespace `Loss`

### Lower-level operators

Advanced users building specialised pipelines can bypass `Min` entirely and call the core operators directly.

#### Trust-region minimisation

`M` is a lower version of `Min`. It is a dyadic operator from which a dyadic function is derived. Usage:

    R←X f M g Y

where `f` is a monadic evaluation function, `g` is a monadic function which takes as argument an `iter cost rel dnorm p` vector and gets called before every convergence check, `Y` is either a two elements vector with the initial guess of parameters and normalised damping factor; or a five elements vector with parameters, damping, lower bounds, upper bounds, and margin. `X` is a vector with the configuration parameters `toli tolc tolr tolg dini dinc ddec dmax dmin`.

The return value of `f` will determine the algorithm to use (LM for residual vectors, BFGS for scalars, Newton if Hessian and gradient are provided). If Jacobian matrices or gradient vectors are not returned by `f`, they will be numerically estimated.

The return value `R` is an `iter cost rel dnorm p` vector.

#### Trust-region Newton solver

`Newton` is the core trust-region Newton solver.

    R←X f Newton g Y

It is analogous to the `M` operator but its evaluation function `f` must be an ambivalent function that returns the value of the objective function, Hessian matrix, and gradient. Optionally, it might also return updated parameters at the front. The first iteration, `f` will be called monadically. In successive iterations, the last accepted parameters, cost, Hessian and gradient will be given as left argument.

#### Hessian approximation

The monadic operators `LM` and `BFGS` take an evaluation function `f` and return ambivalent functions to use as left operand for `Newton`.

    c h g←{X}f LM Y
    c h g←{X}f BFGS Y

Composing these operators as `f LM Newton` or `f BFGS Newton` gives a complete solver.

For `LM`, `f` returns either the residuals, the residuals and the Jacobian, or the residuals, Jacobian, loss values and weights. If no Jacobian is provided, it is estimated numerically. If no loss values and weights are provided, squared residuals are used. The derived function returns the cost, the Gauss-Newton matrix $J^T W J$ as the Hessian model, and the gradient $J^T W y$.

For `BFGS`, `f` returns the scalar objective. The gradient is estimated numerically. The derived function returns the cost, the BFGS approximation of the Hessian (initialised to the identity and updated each iteration from successive parameter and gradient differences), and the gradient. Updates that would violate the curvature condition $s_k^T y_k > 0$ are skipped.

#### Bounded box constraints

`Bounded` is a dyadic operator that wraps a Hessian-providing evaluator to enforce box constraints.

    p c h g←{X}(f Bounded g)Y

The left operand `f` is analogous to the same operand in `Newton`. The right operand is a three-element vector with the lower bounds, upper bounds, and interior margin. Lower and upper bounds may contain infinite entries to indicate no bound on that side. If the margin is omitted, it defaults to `⎕CT`.

## Examples

See [tests](APLSource/Test.aplf).
