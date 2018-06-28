import numpy as np
from matplotlib import pyplot as plt

# Constants
length = 50
v_max = 15
delta_t = 0.01
r_cut = 10
particle_count = 9

# Initialise simulation
positions = np.random.uniform(0, length, (particle_count, 2))
velocities = np.random.uniform(-v_max, v_max, (particle_count, 2))
forces = np.zeros((particle_count, 2), dtype=np.float64)

# Function to calculate the forces using the Lennard-Jones Potential
def forces_calc():

    # Create square matrices of x,y coords for particles for i > j
    pos_tiled_x = np.tile(positions[:, 0], (particle_count, 1))
    pos_tiled_y = np.tile(positions[:, 1], (particle_count, 1))

    r_ij = np.empty((particle_count, particle_count, 2), dtype=np.float64)
    # Finding the position vectors r_ij
    r_ij[:, :, 0] = np.subtract(np.triu(np.transpose(pos_tiled_x)), np.triu(pos_tiled_x))
    r_ij[:, :, 1] = np.subtract(np.triu(np.transpose(pos_tiled_y)), np.triu(pos_tiled_y))
    norm_r_ij = np.linalg.norm(r_ij, axis=2)  # Finding the norm of the position vectors

    index_array = np.less(norm_r_ij, r_cut)  # Indexing array for particles satisfying the threshold condition
    calc_forces = np.zeros((particle_count, particle_count), dtype=np.float64)
    calc_forces[index_array] = np.subtract(np.multiply(np.power(norm_r_ij[index_array], -14), 48),
                                           np.multiply(np.power(norm_r_ij[index_array], -8), 24))  # Potential for ij
    calc_forces[np.isnan(calc_forces)] = 0.0  # Handling NaNs in the resulting upper triangular matrix
    r_ij[:, :, 0] = np.multiply(r_ij[:, :, 0], calc_forces)  # Resulting force in the x component for ij
    r_ij[:, :, 1] = np.multiply(r_ij[:, :, 1], calc_forces)  # Resulting force in the y component for ij
    r_ij[:, :, 0] = np.add(r_ij[:, :, 0], -np.transpose(r_ij[:, :, 0]))  # Corresponding lower triangular entry for x
    r_ij[:, :, 1] = np.add(r_ij[:, :, 1], -np.transpose(r_ij[:, :, 1]))  # Corresponding lower triangular entry for y

    global forces
    forces = np.sum(r_ij, axis=1)  # Resultant net force on each particle i

# Function to evaluate the motion of the particles
def velocity_verlet():

    global positions, velocities
    positions = np.add(np.add(positions, np.multiply(delta_t, velocities)), np.multiply(0.5*(delta_t**2), forces))  # Update position of particles

    # Boundary reflection conditions
    index_array_L = np.greater(positions, length)
    positions[index_array_L] = np.add(-positions[:, 0:2][index_array_L], 2*length)  # For beyond length L
    index_array_0 = np.less(positions, 0)
    positions[index_array_0] = -positions[index_array_0]  # For less than 0

    forces_prev = np.copy(forces)  # Copy forces from previous iteration
    forces_calc()  # Calculate the new forces

    velocities = np.add(velocities, np.multiply(delta_t/2, np.add(forces_prev, forces)))  # Update velocities

    # Boundary velocity inversion conditions
    velocities[np.logical_or(index_array_L, index_array_0)] = -velocities[np.logical_or(index_array_L, index_array_0)]

# Main code
plt.ion()
plt.show()

forces_calc()  # Calculating the initial forces to start the simulation

while True:
    velocity_verlet()
    plt.xlim((0, length))
    plt.ylim((0, length))
    plt.autoscale(False)
    plt.scatter(positions[:, 0], positions[:, 1])  # Plotting
    plt.pause(0.001)
    plt.clf()
