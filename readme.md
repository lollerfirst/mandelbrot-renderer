# Mendelbrot Rendered

This program creates a blank canvas and draws the points which belong to the _MendelBrot Set_ $M$.
Specifically, a point $\mathbf{P}: (x, y)$ translates to a complex value of $c = x + iy$, without considering the scaling.
A point belongs to the Mendelbrot set if the series 
```math
  z(n) = \begin{cases} 0, & \text{if } n-1 = 0\\z(n-1)^2 + c, &\text{if } n-1 > 0\end{cases}
```
converges.
That is $$ c \belong \mathbf{M} \iff lim_{x \to \infty}{z(x)} = k}, & k \belong \mathbb{C} $$ 
