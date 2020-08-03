---
layout: post
title:  "Research on the solution of a type of SDE"
categories: math research
tags: math SDE Fokker–Planck-equation proof
author: Ziyu Zhong
mathjax: true
---

* content
{:toc}

## Introduction

The problem comes from a stochastic PDE:

$$\dfrac{du}{dt}=\Delta u-Eu^2+\xi$$

Where &nbsp; $u$ &nbsp; is a function of $t$ and $x$, $x\in$ [0,1], and &nbsp; $\xi$ &nbsp; is space-time white noise.

We want to understand invariant measure of the stochastic PDE.

However, since it's complicated to deal with SPDE, we can try to make a finite difference approximation. Let’s take n points from [0,1].

Then it becomes an SDE system:

$$dX_t^i=(X_t^{i-1}+X_t^{i+1}-2X_t^{i}-E{X_t^i}^2\cdot X_t^i)dt+dB_t^i$$

where $i=1,...,n$. We take a periodic boundary condition here, namely $X_t^0=X_t^n$ and $X_t^{n+1}=X_t^1$.

We want to figure out the invariant measure of this SDE system.


## Preliminary

### Fokker–Planck equation

The Fokker–Planck equation is a partial differential equation that describes the time evolution of the probability density function of the random variebles. For a general case, if

$$dX_t=\mu (X_t, t)dt+\sigma (X_t, t) dB_t $$

where $X_t$ and $\mu(X_t,t)$ are $n$-dimensional random vectors, $\sigma (X_t,t)$ is an $n\times m$ matrix and $B_t$ is an $m$-dimensional standard Wiener process, the probability density $p(x,t)$ for $X_t$ satisfies the Fokker–Planck-equation:

$$\frac{\partial p(\mathbf{x}, t)}{\partial t}=-\sum\limits_{i=1}^{N} \frac{\partial}{\partial x_{i}}\left[\mu_{i}(\mathbf{x}, t) p(\mathbf{x}, t)\right]+\sum_{i=1}^{N} \sum_{j=1}^{N} \frac{\partial^{2}}{\partial x_{i} \partial x_{j}}\left[D_{i j}(\mathbf{x}, t) p(\mathbf{x}, t)\right]$$

with drift vector $\mu = (\mu_1,...,\mu_N)$ and diffusion tensor $D = \frac{1}{2}\sigma \sigma^T$, i.e.

$$D_{i j}(\mathbf{x}, t)=\frac{1}{2} \sum_{k=1}^{M} \sigma_{i k}(\mathbf{x}, t) \sigma_{j k}(\mathbf{x}, t)$$

For our problem, 

$\mu =$
$\begin{pmatrix}
-2-E{X_t^1}^2 & 1 & 0 & \cdots & 0 & 1 \\\ 
1 & -2-E{X_t^2}^2 & 1 & \cdots & 0 & 0 \\\ 
0 & 1 & -2-E{X_t^3}^2 & 1 & \cdots & 0 \\\ 
\vdots & \vdots & \ddots & \ddots & \ddots &\vdots \\\ 
0 & \cdots & 0 & 1 & -2-E{X_t^{n-1}}^2 & 1 \\\ 
1 & 0 & \cdots & 0 & 1 & -2-E{X_t^n}^2 \\\ 
\end{pmatrix}$

&nbsp;

$\sigma = I$ &nbsp; (identity matrix)

&nbsp;

If our $\mu$ is a constant, the SDE would be an Ornstein–Uhlenbeck process, whose invariant measure is very clear for us. (one can see the derivation from <https://doi.org/10.1186/s13662-019-2214-1>)

However, our $\mu$ is relevant with $t$, since $E{X_t^i}^2$ is a function of $t$. Thus it's tricky to deal with it.

## Try to solve it

....