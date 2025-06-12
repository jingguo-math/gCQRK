# Runge--Kutta generalized Convolution Quadrature for sectorial problems

Code accompanying the paper "Runge--Kutta generalized Convolution Quadrature for sectorial problems" by Jing Guo, Maria Lopez-Fernandez.

## Abstract

We study the application of the generalized convolution quadrature (gCQ) based on Runge--Kutta methods to approximate the solution of an important class of sectorial problems. The gCQ generalizes Lubich's original convolution quadrature (CQ) to variable steps. High-order versions of the gCQ have been developed in the last decade, relying on certain Runge--Kutta methods. The Runge--Kutta-based gCQ has been studied so far in a rather general setting, which includes applications to boundary integral formulations of wave problems. The available stability and convergence results for these new methods are suboptimal compared to those known for the uniform-step CQ, both in terms of convergence order and regularity requirements of the data. Here we focus on a special class of sectorial problems and prove that in these important applications it is possible to achieve the same order of convergence as for the original CQ, under the same regularity hypotheses on the data, and for very general time meshes. In the particular case of data with some known algebraic type of singularity, we also show how to choose an optimally graded time mesh to achieve convergence with maximal order, overcoming the well-known order reduction of the original CQ in these situations. 
An important advantage of the gCQ method is that it allows for a fast and memory-efficient implementation. We describe how the fast and oblivious Runge--Kutta-based gCQ can be implemented and illustrate our theoretical results with several numerical experiments.

## Citation
