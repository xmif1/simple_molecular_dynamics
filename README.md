# simple_molecular_dynamics

Suppose we want to simulate the motion and dynamic interactions between *n* particles in a closed 2-dimensional container. The task at hand is to numerically integrate Newton’s equations for a set of particles, with the interaction force between any two particles being derived from the Lennard Jones potential. The net force on each particle is then to be used to evaluate the system over time, using the Velocity Verlet algorithm for a timestep Delta *t*.

The program was essentially implemented using two functions, one to calculate the resultant force on each particle and the other to carry out the Velocity Verlet algorithm. Heavy use of Numpy arrays and advanced indexing was made, effectively resulting in an implementation without the
use of any looping constructs (besides for plotting). In order to reduce the number of calculations needed to find all of the forces, an upper triangular matrix of relative position vectors was constructed, and the force was then calculated on the position vectors within this upper triangular matrix.

Using Newton’s Third Law of Motion, the corresponding lower triangular matrix could be found (by using the results in the upper triangular matrix rather than carrying out another set of calculations). Summing these two matrices resulted in a matrix with the component forces for the *i* th particle in the *i* th row. Thus by summing over the 1st axis, the resultant force for the *i* th particle was found, where *i* is in {1, ..., *n*}. To further reduce the number of calculations, a norm cutoff parameter was added.

A number of readily available Numpy functions were used for efficiency, including: 'np.transpose' for transposing arrays, 'np.linalg.norm' for finding the position vectors, 'np.random.uniform' for generating the particle positions and velocities, as well as the Numpy functions for array addition and multiplication. Indexing arrays were also generated using 'np.greater' and 'np.less'.

The same holds for the Velocity Velocity algorithm - the required equations were essentially implemented as the addition of arrays, eliminating the need for any 'for' loops. Within the Velocity Verlet algorithm a set of boundary reflection conditions were implemented for particles crossing the container boundary - effectively the particles were reflected off of the boundary, with the assumption of conservation of energy.

Lastly, for the graphical aspect of the program the scatter functionality provided by 'matplotlib' was used to plot the position vectors for each particle with each iteration of the 'while' loop.
