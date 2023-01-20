# Solving the Thomson problem using good ol' observation.
# Of course atoms don't work like that but we can model
# an imaginary system that does work like that.
# Still, I will write the comments as if we are trying to
# solve for electron positions on a unit sphere in 3D space.
import random
import matplotlib.pyplot as plt
import numpy as np

from vector3 import *

class electron:
    def __init__(self, pos):
        self.pos = pos
        self.vel = vec3(0,0,0)

    def calc_accel(self, other_electron):
        accel_dir = (self.pos - other_electron.pos).normalized()

        # Here, we won't bother calculating an actual coulomb force.
        # The magnitude of the force doesn't matter if we are merely
        # trying to find the state at which the potential energy is
        # minimum. So the numerator is just set to 1 to reduce
        # computational expense.
        accel_mag = 1 / ((self.pos - other_electron.pos).mag()**2)
        return accel_dir * accel_mag

    def update_vel(self, accel):
        self.vel += accel

    def update_pos(self):
        self.pos += self.vel

# takes:
# n = number of electrons
# tolerance = position convergence tolerance
def solve_thomson(n, tolerance, relaxation_factor=1):
    # Initialize electrons at random positions on a unit sphere
    # that is centered at position (0,0,0).
    electrons = []
    for i in range(n):
        new_pos = vec3()
        new_pos.x = random.uniform(-1, 1)
        new_pos.y = random.uniform(-1, 1)
        new_pos.z = random.uniform(-1, 1)
        new_pos = new_pos.normalized()
        
        new_electron = electron(new_pos)
        electrons.append(new_electron)

    # Electrons are initialized. Now we will let them push
    # each other electromagnetically.
    converged = False
    cycle_num = 0
    while not converged:
        current_electron_accels = []
        # we will calculate acceleration of each electron...
        for e in electrons:
            total_electron_accel = vec3(0, 0, 0)

            # ...due to each electron...
            for e2 in electrons:
                
                # ...that is not the electron we are
                # currently calculating the acceleration of.
                if not e == e2:
                    total_electron_accel += e.calc_accel(e2)

            current_electron_accels.append(total_electron_accel)

        # Accelerations have been calculated. Now, update their velocities.
        for idx_e in range(len(electrons)):
            e = electrons[idx_e]
            # apply 10 percent "drag" so that electrons lose velocity over steps and do not bounce around their lowest energy points
            e.update_vel((current_electron_accels[idx_e] - current_electron_accels[idx_e] * 0.9) * relaxation_factor)

            # limit velocity components so that the electrons only move on
            # the surface of the sphere
            radial_vel = e.vel.dot(e.pos.normalized())
            e.vel = e.vel - e.pos.normalized() * radial_vel

        # Velocities are updated. Now, update positions.
        old_positions = []
        new_positions = []
        for e in electrons:
            old_positions.append(e.pos)
            e.update_pos()
            e.pos = e.pos.normalized() # limit movement domain to sphere surface
            new_positions.append(e.pos)

        # check for position convergence
        convergence_list = []
        for idx_pos in range(len(old_positions)):
            old_pos = old_positions[idx_pos]
            new_pos = new_positions[idx_pos]

            convergence_list.append((new_pos - old_pos).mag())

        converged = all(i < tolerance for i in convergence_list)
        cycle_num += 1

        if cycle_num % 2500 == 0:
            print("Convergence:", convergence_list)

    return new_positions

n = 13 # number of electrons
tolerance = 1e-2
relaxation_factor = 1 # try playing with this value if it takes long to converge
final_positions = solve_thomson(n, tolerance, relaxation_factor)

# print solution numerically
print("\nSolution:")
for electron_pos in final_positions:
    print(electron_pos)

# plot solution
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

# first, draw a unit sphere
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
i = 0.9 * np.outer(np.cos(u), np.sin(v))
j = 0.9 * np.outer(np.sin(u), np.sin(v))
k = 0.9 * np.outer(np.ones(np.size(u)), np.cos(v))
ax.plot_surface(i, j, k)

# re-format solution into some form pyplot can utilize
xs = []
ys = []
zs = []
for final_pos in final_positions:
    xs.append(final_pos.x)
    ys.append(final_pos.y)
    zs.append(final_pos.z)

# plot electrons
ax.scatter(xs, ys, zs, color="r")
plt.title("Solution of Thomson problem for n=" + str(n) + " electrons")

plt.show()
