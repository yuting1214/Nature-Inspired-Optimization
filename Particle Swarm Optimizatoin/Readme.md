1. Particle Swarm Optimization
  + Introduction
  + Pseudo-code
  + Configuration
  + Experimental Result
  
  \newpage

# 1. Particle Swarm Optimization
## (1-1) Introduction
Particle swarm optimization, or PSO, was developed by Kennedy and Eberhart in 1995 $[1]$ and has become one of the most widely used swarm-intelligence-based algorithms due to its simplicity and flexibility.

PSO has roots in two main component methodologies. Perhaps more obvious are its ties to artificial life in general, and to bird flocking, fish schooling, and swarming theory in particular. It is also related, however, to evolutionary computation, and has ties to both genetic
algorithms and evolutionary programming.

The optimization's metaphor had two cognitive aspects, individual learning and learning from a social group.
Where an individual finds itself in a problem(or search) space it can use its own experience and that of its peers to move itself toward the solution.
$$ v_{i}^{t+1} = v_{i}^{t}  +\alpha \varepsilon_{1}(p_{i} - x_{i})+\beta \varepsilon_{2}(p_{g} - x_{i}) \tag{1}$$
$$  x_{t+1} = x_{t} + v_{t+1} \tag{2}$$
,where constants $\alpha$ and $\beta$ determine the balance between the influence of the individual knowledge $\alpha$ and that of the group $\beta$ (both set initially to 2), $\varepsilon_{1}$ and $\varepsilon_{2}$ are uniformly distributed random numbers defined by some upper limit, $\varepsilon_{max}$,that is a parameter of the algorithm, $p_{i}$ and $p_{g}$ are respectively the individual's current best position and the global current best position, and $x_{i}$ is the current position in the dimension considered.

However, by far the most problematic characteristic of PSO is its propensity to converge, prematurely, on early best solutions. Many strategies have been developed in attempts to overcome this situation. Among them, the most popular method are inertia and constriction. The inertia term, $\omega$, was introduced thus (Shi and Eberhart 1998):
$$ v_{i}^{t+1} = \omega v_{i}^{t}  +\alpha \varepsilon_{1}(p_{i} - x_{i})+\beta \varepsilon_{2}(p_{g} - x_{i}) \tag{3}$$

