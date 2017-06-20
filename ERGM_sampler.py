import numpy as np

def ERGM_sampler(dg_funs,thetas,A0,K=100,burn_in=0.1):
    """Generate a list of samples from an ERGM model, using Gibbs Sampling.
    
    The sampler iterates over all possible edges in the graph. If G is the 
    graph as it is, and G' is the graph with the current edge flipped, then 
    the edge is flipped with probability given by the Hasting ratio
    
        H = min(1,P(G')/P(G))
        
    In an ERGM, the quantity P(G')/P(G) is the dot product of theta (a vector
    of parameters) and dg (a vector of changes in graph metrics).
    
    Parameters
    ----------
    dg_funs : list of functions
        List of functions that calculate the change in graph metrics due to
        the flipping of a given edge. Must take as inputs (A,i,j), i.e. the
        full graph, as well as the particular edge to be investigated. Must
        return a float.
        
    thetas : array
        Array of parameters, telling us the relative importance of the dg
        functions.
        
    A0 : array
        Initial array. Assumed to be a simple graph, in that it is both 
        unweighted and undirected, with no self-loops.
        
    K : int, optional (default=100)
        Size of ensemble. 
        
    burn_in : float, optional (default=0.1)
        Fraction of the ensemble to be discarded after sampling.
        
    Returns
    -------
    A_list : list of arrays
        Samples from the ERGM distribution. Length is ceil(N*(1-burn_in)).
        
    Notes
    -----
    See README.md for more details on the mathematics. This function iterates over
    the upper triangular portion of the array, with the bottom portion automatically
    assigned to be equal to the top, and the diagonal ignored.
    """
    n,m = A0.shape
    # Check inputs for validity
    if n != m:
        raise ValueError('Initial array must be square.')
    if len(dg_funs) != len(thetas):
        raise ValueError('Number of parameters must match number of metrics.')
    # if a non-symmetric matrix is input, we just ignore the lower portion
    A = np.triu(A0) + np.triu(A0).T
    A_list = []
    burn_in_len = int(K*burn_in) # truncates towards zero
    for k in range(K):
        for i in range(n):
            for j in range(i+1,n):
                dg_vec = np.array([fun(A,i,j) for fun in dg_funs])
                logH = min(0,sum(thetas*dg_vec)) # log of Hastings Ratio
                if np.random.binomial(1,np.exp(logH)):
                    A[i,j] = 1-A[i,j] # flips the edge
                    A[j,i] = A[i,j] # ensure symmetry
        if k >= burn_in_len:
            A_list.append(A.copy())
    return A_list
