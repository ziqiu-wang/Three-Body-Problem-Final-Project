import numpy as np
from matplotlib import pyplot as plt
from projectIntegrator import Body, solve_n_body_problem, plot_helper, plot_solution

# Testing different orbit shapes
bodies = [Body(2, np.array([0.0,0.0,0.0]), np.array([-1.0,0.0,0.0])),
          Body(2, np.array([0.0,1.0,0.0]), np.array([1.0,0.0,0.0]))]
plot_solution(bodies, 0.0001, 50000) # circular orbits

bodies = [Body(2, np.array([0.0,0.0,0.0]), np.array([-1.0,0.0,0.0])),
          Body(2, np.array([0.0,1.4,0.0]), np.array([1.0,0.0,0.0]))]
plot_solution(bodies, 0.0001, 150000) # elliptical orbits

bodies = [Body(2, np.array([0.0,0.0,0.0]), np.array([-1.0,0.0,0.0])),
          Body(5, np.array([0.0,1.2,0.0]), np.array([1.0,0.0,0.0]))]
plot_solution(bodies, 0.0001, 50000) # elliptical orbits

bodies = [Body(1, np.array([0.0,1.0,0.0]), np.array([1.0,0.0,0.0])),
          Body(7, np.array([0.0,0.0,0.0]), np.array([0.0,0.0,0.0]))]
plot_solution(bodies, 0.0001, 50000) # elliptical orbits

bodies = [Body(1, np.array([0.0,0.0,0.0]), np.array([1.0,0.0,0.0])),
          Body(2, np.array([1.0,1.0,0.0]), np.array([-1.0,0.0,0.0]))]
plot_solution(bodies, 0.0005, 400000) # elliptical orbits, slanted


# Testing periods of two body systems
bodies = [Body(2, np.array([0.0,0.0,0.0]), np.array([-0.03,0.0,0.0])),
          Body(0.002, np.array([0.0,1.0,0.0]), np.array([1.0,0.0,0.0]))]
solve_n_body_problem(bodies, 0.0001, 54000)

bodies = [Body(2, np.array([0.0,0.0,0.0]), np.array([-1.0,0.0,0.0])),
          Body(2, np.array([0.0,1.0,0.0]), np.array([1.0,0.0,0.0]))]
solve_n_body_problem(bodies, 0.0001, 200000)


# Testing convergence of time step
dt_lst = [0.5, 0.2, 0.1, 0.01, 0.001, 0.0001]
theta = np.linspace(0, 2 * np.pi, 100000)
x_0 = 0.5 * np.cos(theta)
y_0 = 0.5 * np.sin(theta)
for dt in dt_lst:
    plt.plot(x_0, y_0, "olive")
    bodies = [Body(2, np.array([0.0, 0.0, 0.0]), np.array([-1.0, 0.0, 0.0])),
              Body(2, np.array([0.0, 1.0, 0.0]), np.array([1.0, 0.0, 0.0]))]  # period: pi
    flat_r_arr = solve_n_body_problem(bodies, dt, int(50.0/dt))  # use this, try dt = 0.1, 0.01, 0.001, 0.0001, (0.00005)
    plot_helper(bodies, flat_r_arr)
    plt.title(f"$dt = {dt}$")
    plt.show()