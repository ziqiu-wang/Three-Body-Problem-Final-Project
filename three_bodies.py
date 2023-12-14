## This file is used to produce the plots for three-bodies systems, in the CM frame.


import numpy as np
import matplotlib.pyplot as plt
from projectIntegrator import solve_n_body_problem, Body, plot_helper, plot_solution


## PART 3 Q1

bodies = [Body(2, np.array([0.0,0.0,0.0]), np.array([-1.0,0.0,0.0])),
          Body(2, np.array([0.0,1.0,0.0]), np.array([1.0,0.0,0.0])),
          Body(0.01, np.array([1.0,-1.0,0.0]), np.array([-1.0,1.0,0.0]))]
plot_solution(bodies, 0.0001, 80000)

bodies = [Body(2, np.array([0.0,0.0,0.0]), np.array([-1.0,0.0,0.0])),
          Body(2, np.array([0.0,1.0,0.0]), np.array([1.0,0.0,0.0])),
          Body(0.01, np.array([-0.4,-1.2,0.0]), np.array([-2.0,0.8,0.0]))]
plot_solution(bodies, 0.0001, 80000)    # another set of initial conditions

## PART 3 Q4

bodies = [Body(2, np.array([0.0,0.0,0.0]), np.array([0.0,0.0,0.0])),
          Body(0.02, np.array([0.0,2.0,0.0]), np.array([1.0,0.0,0.0]))]
plot_solution(bodies, 0.0001, 200000)   # just the two particles

r_lst = [0.25, 0.5, 1.0, 1.7, 1.9, 1.97, 1.99, 1.997, 1.999,
         2.0, 2.001, 2.004, 2.01, 2.03, 2.1, 2.2, 2.5, 3.0, 5.0, 10.0]
for r in r_lst:
    bodies = [Body(2, np.array([0.0,0.0,0.0]), np.array([0.0,0.0,0.0])),
              Body(0.02, np.array([0.0,2.0,0.0]), np.array([1.0,0.0,0.0])),
              Body(0.02, np.array([r,0.0,0.0]), np.array([0.0,-np.sqrt(2/r),0.0]))]
    flat_r_arr = solve_n_body_problem(bodies, 0.001, 100000)
    plot_helper(bodies, flat_r_arr)
    plt.title(f"$r = {r}$")
    plt.show()

##L1 Lagrange Points
bodies = [Body(50, np.array([-1.0,0.0,0.0]), np.array([0.0,3.8,0.0])),
          Body(50, np.array([1.0,0.0,0.0]), np.array([0.0,-3.8,0.0])),
          Body(2, np.array([0.0,0.0,0.0]), np.array([0.0,0.0,0.0]))]

flat_r_arr = solve_n_body_problem(bodies, 0.0001, 50000)

plot_helper(bodies, flat_r_arr)
plt.title('L1 in Rotating Frame')
plt.show()

##L2/L3 Lagrange Points
bodies = [Body(2, np.array([0.0,-0.5,0.0]), np.array([-1.0,0,0.0])),
          Body(2, np.array([0.0,0.5,0.0]), np.array([1.0,0.0,0.0])),
          Body(.01, np.array([0.0,-.5992,0.0]), np.array([1.1984,0.0,0.0]))]
         

flat_r_arr = solve_n_body_problem(bodies, 0.0001, 50000)

plot_helper(bodies, flat_r_arr)
plt.title('L2/L3 in Rotating Frame')
plt.show()
   
##L4/5 Lagrange Points
bodies = [Body(25, np.array([0.0,0.0,0.0]), np.array([0.0,0.0,0.0])),
          Body(1, np.array([1.0,0.0,0.0]), np.array([0.0,5.099,0.0])),
          Body(.01, np.array([-0.4615,0.866,0.0]), np.array([-4.33,-2.36,0.0]))
          ]



flat_r_arr = solve_n_body_problem(bodies, 0.0001, 50000)

plot_helper(bodies, flat_r_arr)
plt.title('L4/L5 in Rotating Frame')
plt.show()


## Figure eight

# Define the masses and initial conditions for the figure-eight solution
mass = 1  # for equal masses
# Positions
pos1 = np.array([-0.97000436, 0.24308753, 0])
pos2 = np.array([0, 0, 0])
pos3 = np.array([0.97000436, -0.24308753, 0])

vel1 = np.array([0.466203685, 0.43236573, 0])
vel2 = np.array([-0.93240737, -0.86473146, 0])
vel3 = np.array([0.466203685, 0.43236573, 0])

# Create the bodies
body1 = Body(mass, pos1, vel1)
body2 = Body(mass, pos2, vel2)
body3 = Body(mass, pos3, vel3)

# List of bodies
bodies = [body1, body2, body3]


    
    
flat_r_arr = solve_n_body_problem(bodies, 0.0001, 50000)

plot_helper(bodies, flat_r_arr)
plt.title('Three Bodies of Similar Masses')
plt.show()
