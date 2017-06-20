# Gibbs Sampler for Exponential Random Graph Models

### Peter Wills, 6/20/2017

## Mathematical Background

The following discussion assumes a basic understanding of network science. Please see [1] for the relevant background.

This module consists solely of a Gibbs sampler that can generate draws from an Exponential Random Graph Model (ERGM). This model assigns each graph G  a probability

$$ P(G) = \frac{1}{Z} \exp\left\{\sum_{i=1}^r \theta_i g_i(G)\right\}$$

where $g_i(G)$ is some graph metric (e.g., $g_1(G)$ is the volume of $G$ while $g_2(G)$ is the number of triangles) and the $\theta_i$ are parameters that tell us "how important" each metric is in the distribution.

The trick here is that we don't want to, and in practice cannot, calculate the normalizing factor $Z$. So we need a sampling method that doesn't depend on $Z$. This is where Gibbs sampling comes in.

A detailed discussion of this method of sampling  outline the algorithm. We'll denote subsequent graphs with $G^t$, and edges in a graph with $G_{ij}$. 

- Start with a graph $G^0$.
- For $t=1,2,\ldots$:
    - For $i=1,\ldots,n$:
        - For $j=i+1,\ldots,n$:
            - Let $G^*$ be the graph where $G^*_{ij} = 1-G^t_{ij}$ and all other edges are the same.
            - Calculate the Hastings Ratio:
            $$ H = \min\left\{1,\frac{P(G^*)}{P(G^t)}\right\} = \min\left\{ 1, \exp\left\{\sum_{i=1}^r \theta_i (g_i(G^*)-g_i(G^t))\right\}\right\}$$
            - With probability $H$, assign $G^{t+1}=G^*$. Otherwise, assign $G^{t+1}=G^t$.
            
This process loops through all the ${n\choose 2}$ edges in the graph, so will be very slow unless the calculation of the terms $\Delta g_i = g_i(G^*)-g_i(G^t)$ is quick. This depends on our choice of metric $g_i$, of course. We want to focus on ones that have easy update formulas, that can ideally be computed in linear time.

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
	