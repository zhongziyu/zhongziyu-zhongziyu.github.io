---
layout: post
title:  "Research on the solution of a type of SDE"
categories: math research
tags: math SDE Fokker–Planck-equation proof
author: Ziyu Zhong
---

* content
{:toc}

# 1 Introduction

The problem comes from a stochastic PDE:

$$\dfrac{du}{dt}=\Delta u-Eu^2+\xi$$

Where &nbsp; $u$ &nbsp; is a function of $t$ and $x$, $x\in$ [0,1], and &nbsp; $\xi$ &nbsp; is space-time white noise.

We want to understand invariant measure of the stochastic PDE.

However, since it's complicated to deal with SPDE, we can try to make a finite difference approximation. Let’s take n points from [0,1].

Then it becomes an SDE system:

$$dX_t^i=(X_t^{i-1}+X_t^{i+1}-2X_t^{i}-E{X_t^i}^2\cdot X_t^i)dt+dB_t^i$$

where $i=1,...,n$. We take a periodic boundary condition here, namely $X_t^0=X_t^n$ and $X_t^{n+1}=X_t^1$.

We want to figure out the invariant measure of this SDE system.


# 2 Preliminary

## 2.1 Fokker–Planck equation

The Fokker–Planck equation is a partial differential equation that describes the time evolution of the probability density function of the random variebles. For a general case, if

$$dX_t=\mu (X_t, t)dt+\sigma (X_t, t) dB_t $$

where $X_t$ and $\mu(X_t,t)$ are $n$-dimensional random vectors, $\sigma (X_t,t)$ is an $n\times m$ matrix and $B_t$ is an $m$-dimensional standard Wiener process, the probability density $p(x,t)$ for $X_t$ satisfies the Fokker–Planck-equation:

$$\frac{\partial p(\mathbf{x}, t)}{\partial t}=-\sum\limits_{i=1}^{N} \frac{\partial}{\partial x_{i}}\left[\mu_{i}(\mathbf{x}, t) p(\mathbf{x}, t)\right]+\sum_{i=1}^{N} \sum_{j=1}^{N} \frac{\partial^{2}}{\partial x_{i} \partial x_{j}}\left[D_{i j}(\mathbf{x}, t) p(\mathbf{x}, t)\right]$$

with drift vector $\mu = (\mu_1,...,\mu_N)$ and diffusion tensor $D = \frac{1}{2}\sigma \sigma^T$, i.e.

$$D_{i j}(\mathbf{x}, t)=\frac{1}{2} \sum_{k=1}^{M} \sigma_{i k}(\mathbf{x}, t) \sigma_{j k}(\mathbf{x}, t)$$

For our problem, 

$\mu =$
$\begin{pmatrix}
-2-E{X_t^{(1)}}^2 & 1 & 0 & \cdots & 0 & 1 \\\ 
1 & -2-E{X_t^{(2)}}^2 & 1 & \cdots & 0 & 0 \\\ 
0 & 1 & -2-E{X_t^{(3)}}^2 & 1 & \cdots & 0 \\\ 
\vdots & \vdots & \ddots & \ddots & \ddots &\vdots \\\ 
0 & \cdots & 0 & 1 & -2-E{X_t^{(n-1)}}^2 & 1 \\\ 
1 & 0 & \cdots & 0 & 1 & -2-E{X_t^{(n)}}^2 \\\ 
\end{pmatrix}$

&nbsp;

$\sigma = I$ &nbsp; (identity matrix)

&nbsp;

If our $\mu$ is a constant, the SDE would be an Ornstein–Uhlenbeck process, whose invariant measure is unique and very clear for us. (one can see the derivation from <https://doi.org/10.1186/s13662-019-2214-1>)

However, our $\mu$ is relevant with $t$, since $E{X_t^{(i)}}^2$ is a function of $t$. Thus it's tricky to deal with it.

# 3 Try to solve it

## 3.1

We assume $X_t$ is stationary, then $E{X_t^{(i)}}^2$ is constant. Let the constant be $u_i$. Then we can regard $\mu$ as a constant. So the SDE becomes an Ornstein–Uhlenbeck process. We can obtain that its invariant measure should be a gaussian distribution $X_t\sim N(0,-\frac{1}{2}\mu^{-1}) $. For completeness, I give a proof as follows.

For our SDE with the initial value $\mathbf{X_0}=\mathbf{x_0}$, the Fokker-Planck equation is:

$\dfrac{\partial p}{\partial t} = -\sum\limits_{i=1}^{n}\dfrac{\partial}{\partial x_i}[(x_{i-1}+x_{i+1}-(2+u_i)x_i)p]+\dfrac{1}{2}\Delta p \qquad (x_0:=x_n,x_{n+1}:=x_1)$

$p(\mathbf{x},0) = \delta (\mathbf{x}-\mathbf{x^0}) \qquad \quad$ ( initial condition )

We take the $n$-dimensional Fourier transform of the equations, we get:

$\dfrac{\partial \hat p}{\partial t} = \sum\limits_{i=1}^{n}\lambda_i[\dfrac{\partial \hat p}{\partial \lambda_{i-1}}+\dfrac{\partial \hat p}{\partial \lambda_{i+1}}-(2+u_i)\dfrac{\partial \hat p}{\partial \lambda_{i}}]-\dfrac{1}{2}\lambda^T \lambda\ \hat p $

$\hat p(\lambda ,0) = exp(-i\lambda^T\mathbf{x^0})$

where $\hat p(\lambda ,t) $ is the $n$-dimensional Fourier transform of $p(\mathbf{x},t)$.

&nbsp;

Note that the equation is a first-order partial differential equation, so we will apply the method of characteristics.

(For the convienience of caculation, we let $A:=-\mu$)

Consider the system:

$\dfrac{d\lambda}{dt} = A\lambda$

with initial condition $\lambda(0) = C$. The solution of this system is:

$\lambda = e^{tA}\cdot C$

Then,

$\dfrac{d\hat p}{dt} = \sum\limits_{i=1}^{n}\dfrac{\partial \hat p}{\partial \lambda_i}\cdot \dfrac{d \lambda_i}{dt} + \dfrac{\partial \hat p}{\partial t}$

$ \quad = -\sum\limits_{i=1}^{n}\lambda_i[\dfrac{\partial \hat p}{\partial \lambda_{i-1}}+\dfrac{\partial \hat p}{\partial \lambda_{i+1}}-(2+u_i)\dfrac{\partial \hat p}{\partial \lambda_{i}}] + \dfrac{\partial \hat p}{\partial t}$

$ \quad =-\dfrac{1}{2}\lambda^T \lambda\ \hat p$

$ \quad =-\dfrac{1}{2}C^Te^{tA^T}e^{tA}C\hat p$

$ \quad =-\dfrac{1}{2}C^Te^{2tA}C\hat p \qquad$ ($A$ is symmetric)

$ \hat p(0) = \hat p(\lambda (0) ,0) = exp(-iC^T\mathbf{x^0})$

Thus,

$\hat p = exp(-iC^T\mathbf{x^0} - \dfrac{1}{2}\displaystyle{ \int_0^t C^Te^{2sA}C ds }  )$

$\quad = exp(-iC^T\mathbf{x^0} - \dfrac{1}{2}C^T(2A)^{-1}[e^{2tA}-I]C)$

$ \overset{C=e^{-tA}\cdot \lambda}{=} exp(-i\lambda^Te^{-tA}\mathbf{x^0} - \dfrac{1}{4}\lambda^Te^{-tA}A^{-1}[e^{2tA}-I]e^{-tA}\lambda)$

$\quad = exp(-i\lambda^Te^{-tA}\mathbf{x^0} - \dfrac{1}{4}\lambda^TA^{-1}[I-e^{-2tA}]\lambda)$

Since the characteristic function is the Fourier transform with opposite sign in the complex exponential,

$X_t\sim N(e^{-tA}\mathbf{x^0}, \dfrac{1}{2}A^{-1}[I-e^{-2tA}])$

And because $A$ is a positive definite matrix when $u_i>0$ (we won't prove it here, since it's tedious), when $t\to \infty$, $e^{-tA}\to \mathbf{0}$, $X_t \sim N(\mathbf{0}, \dfrac{1}{2}A^{-1})$ .

## 3.2

We now know that when $X_t$ is stationary, its distribution would be $N(\mathbf{0}, \dfrac{1}{2}A^{-1})$, 

$A =$
$\begin{pmatrix}
2+u_1 & -1 & 0 & \cdots & 0 & -1 \\\ 
-1 & 2+u_2 & -1 & \cdots & 0 & 0 \\\ 
0 & -1 & 2+u_3 & -1 & \cdots & 0 \\\ 
\vdots & \vdots & \ddots & \ddots & \ddots &\vdots \\\ 
0 & \cdots & 0 & -1 & 2+u_{n-1} & -1 \\\ 
-1 & 0 & \cdots & 0 & -1 & 2+u_n \\\ 
\end{pmatrix}$

$u_i = E{X_t^{(i)}}^2$

However, we need to verify $u_i$'s existence and uniqueness, which is equivalent to verify the existence and uniqueness of the solutions to the equation system:

$$\begin{cases}
u_1=\dfrac{1}{2}A^{-1}_{1,1}\\[2ex]
u_2=\dfrac{1}{2}A^{-1}_{2,2} \\[2ex]
\vdots \\[2ex]
u_n=\dfrac{1}{2}A^{-1}_{n,n} 
\end{cases}$$

($ A^{-1}_{i,j} $ means the $(i,j)$ entry of the matrix $A^{-1}$)

(Note that $\dfrac{1}{2}A^{-1}_{i,i}$ is the variance of the stationary $X_t^{(i)}$)

### 3.2.1 Existence

Let $u:=u_1=u_2=\dots =u_n$, we are going to find a solution of the equations.

$C_n:=$
$\begin{vmatrix}
2+u & -1 & 0 & \cdots & 0 & -1 \\\ 
-1 & 2+u & -1 & \cdots & 0 & 0 \\\ 
0 & -1 & 2+u & -1 & \cdots & 0 \\\ 
\vdots & \vdots & \ddots & \ddots & \ddots &\vdots \\\ 
0 & \cdots & 0 & -1 & 2+u & -1 \\\ 
-1 & 0 & \cdots & 0 & -1 & 2+u \\\ 
\end{vmatrix}_n$

$D_n:=$
$\begin{vmatrix}
2+u & -1 & 0 & \cdots &  0 \\\ 
-1 & 2+u & -1 & \cdots &  0 \\\ 
\vdots & \ddots & \ddots & \ddots  &\vdots \\\ 
0 & \cdots & -1 & 2+u & -1 \\\ 
0 & \cdots & 0 & -1 & 2+u \\\ 
\end{vmatrix}_n$

It's easy to find that $C_n=(2+u)D_{n-1}-2(D_{n-2}+1)$ by expanding $C_n$.

And we can caculate $D_n$ by solving a difference equation, the result is:

$D_n=\dfrac{b^n-a^n}{b-a}$ , $\quad a=\dfrac{2+u-\sqrt{u^2+4u}}{2}$ , $\quad b=\dfrac{2+u+\sqrt{u^2+4u}}{2}$

Since our $u>0$ , then $b>1$ , $a=\dfrac{1}{b}<1$ , $D_n>0$ .

&nbsp;

We see when $u_1=u_2=\dots =u_n=u$ ,

$A_{i,i}^{-1}=\dfrac{D_{n-1}}{C_n}$ ,
since $A^{-1}=\dfrac{A^*}{\mid A\mid}$ .

So the equations become: $u=\dfrac{D_{n-1}}{2C_n}$ .

Both sides multiple $C_n$ and subtitute it by $D_{n-1}$ and $D_{n-2}$, it becomes:

$(2u^2+4u-1)D_{n-1}-4u(D_{n-2}+1)=0$ ,

Let $f(u)$ be the left side of the equation,

$f(u^{1})<0$ , when $u^{1}$ is positive root of $2u^2+4u-1$ , which is $\dfrac{-4+\sqrt{24}}{4}$ .

$f(u^{2})>0$ , when $u^{2}$ is large enough since $D_n$ is $O(u^{n-1})$ .

Thus $f(u)$ has at least one root.


### 3.2.2 Uniqueness

We haven't solved the uniqueness of the solution, but I can show you some ideas that may help us.

####  3.2.2.1 Simple cases

For the simple cases $n=2, 3$ , We can plot a figure of the equations by the computer and find that the equations have only one positive solution:

$n=2$ :

<img src="/imgs/2d.png" width = "500"  alt="2d" align=center />

&nbsp;

$n=3$ :

<img src="/imgs/3d.png" width = "600"  alt="2d" align=center />


