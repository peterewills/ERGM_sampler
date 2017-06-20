# Gibbs Sampler for Exponential Random Graph Models

### Peter Wills

## Mathematical Background

The following discussion assumes a basic understanding of network science. Please see [1] for the relevant background.

This module consists solely of a Gibbs sampler that can generate draws from an Exponential Random Graph Model (ERGM). This model assigns each graph G  a probability

	P(G) = Z^(-1) exp[sum_i(theta_i g_i(G))]

where the `g_i(G)` are some graph metrics (e.g., `g_1(G)` is the volume of G while `g_2(G)` is the number of triangles) and the `theta_i` are parameters that tell us "how important" each metric is in the distribution.

The trick here is that we don't want to, and in practice cannot, calculate the normalizing factor Z. So we need a sampling method that doesn't depend on Z. This is where Gibbs sampling comes in.

We will now outline the algorithm. We use the adjacency matrix `A` as our representation of the graph.
            
	A = A0
	for k = 1,2,...: # iterate as long as you wish
		for i = 1,2,...,n: # size of the matrix
			for j = i+1,i+2,...,n:
				H = min( 1,exp(sum(theta*dg(A,i,j))) )
				with probability H:
					A[i,j] = 1 - A[i,j]
		record A to sample list
            
The array `dg(A,i,j)` is the change in the metrics `g` if we flip edge `i,j`, starting with the state `A`. This process loops through all the `binomial(n,2)` edges in the graph and calculates the change in metrics $dg$ at each step, so will be very slow unless the calculation of the terms `dg(A,i,j)` is quick. This depends on our choice of metrics, of course. We want to focus on ones that have easy update formulas, that can ideally be computed in linear time.

Note that most practicioners recommend using some kind of "burn-in period," where the graphs are not recorded. This is based on the assumption that our algorithm is not sufficiently near convergence at this point. In this module, the burn-in fraction is a user-defined parameter, set at 10% by default.

A more detailed discussion of this sampling method, along with a lot of good information about ERGMs, can be found in [2].

## Using `ERGM_sampler`

Due to the simplicity of the module, no installation procedure has been established; simply place the file `ERGM_sampler.py` in your working directory, or at some location in your `PYTHONPATH`, and then run

	from ERGM_sampler import ERGM_sampler
	
to import the relevant function.

For an example of usage, see `example.ipynb`.

## References

[1] Mark Newman. *Networks: An Introduction.* Oxford University Press, 2010.

[2] Dean Lusher, Johan Koskinen, and Garry Robins. *Exponential Random Graph Models for Social Networks: Theory, Method, and Applications.* Cambridge University Press, 2013
	