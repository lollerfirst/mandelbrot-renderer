# Mendelbrot Rendered

This program creates a blank canvas and draws the points which belong to the _MendelBrot Set_ $M$.
Specifically, a point $\vb{P}: (x, y)$ translates to a complex value of $c = x + iy$, without considering the scaling.
A point belongs to the Mendelbrot set if the series $ z(n) = \begin{cases} 0, & \text{if } n-1 = 0\\z_{n-1}^2 + c, &\text{if } n-1 > 0\end{cases} $
$$ c \belong \vb{M} \iff lim_{x \to \infty}{} $$ 
