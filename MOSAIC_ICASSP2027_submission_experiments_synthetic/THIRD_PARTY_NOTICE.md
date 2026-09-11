# Third-party code notice

The `solver/pgl2021` directories contain MATLAB functions derived from the public implementation accompanying *Product Graph Learning From Multi-Domain Data With Sparsity and Rank Constraints*. The supplied repository snapshot did not contain a license file. Confirm the upstream redistribution terms before publishing these files in a public repository.

The callable `Learn_PGL.m` and `PGL_solver.m` functions expose the paper implementation as ordinary MATLAB functions. The iterative routine returns when the stated residual tolerance is reached; continuing the public script's loop after this point recomputes the same residual without changing the iterate.

The `solver/zw2024` directories contain the implementation used for *Learning Multiplex Graph With Inter-Layer Coupling*. The experiment selects the exact-projection backend for its convex graph subproblem.


