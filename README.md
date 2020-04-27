# A Collection of Planar Contact Models

6 planar rigid-body contact models used in the following papers:

[Empirical evaluation of common contact models for planar impact](https://ieeexplore.ieee.org/abstract/document/7989389)

[Fundamental Limitations in Performance and Interpretability of Common Planar Rigid-Body Contact Models](https://link.springer.com/chapter/10.1007/978-3-030-28619-4_41)

Implementation papers are documented in the model source codes.

## Usage:

```
Inputs:
M := denotes the  3x3 inertia matrix
s := the contact tangent vector (e.g. [1, 0, rx*sin(theta)]) a.k.a. the tangent element of the contact Jacobian
v := velocity of the center of mass in the world frame
ha := sim time step
mu := coefficient of friction
ep := coefficient of restitution
```
```
Outputs:
v_plus := Next time-step velocity
z := Computed impulses (normal and the two tangent in the plane)
```

