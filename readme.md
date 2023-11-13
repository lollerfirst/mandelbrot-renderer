# Mendelbrot Renderer

This program creates a blank canvas and draws the points which belong to the _MendelBrot Set_ $\mathbf{M}$.

Specifically, a point $\mathbf{P}: (x, y)$ translates to a complex value of $c = x + iy$, without considering the scaling.

A point belongs to the Mendelbrot set if the series 
```math
  z_c(n) := \begin{cases} 0, & \text{if } n-1 = 0\\z_c(n-1)^2 + c, &\text{if } n-1 > 0\end{cases}
```
is bounded.
That is if:
```math
 c \in \mathbf{M} \iff \mathbb{R}(lim_{x \to \infty}{z_c(x)}) \leq \mathbb{R}(k) \land \mathbb{I}(lim_{x \to \infty}{z_c(x) k})  \leq  k\ \  k \in \mathbb{C}
``` 
