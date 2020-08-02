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

$$dX_i=(X_{i-1}+X_{i+1}-2X_{i}-EX_i^2\cdot X_i)dt+dB_i$$

where i=1,...,n. We take a periodic boundary condition here, namely $X_0=X_n$ and $X_{n+1}=X_1$.

We want to figure out the invariant measure of this SDE system.


## Fokker–Planck equation

....

## Try to solve it

....