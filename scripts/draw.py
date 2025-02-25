"""
=======================================================

Python script for instance solution visualization
Input/Output should be changed manually for every instance like the verifier

=======================================================
"""


import json
import matplotlib.pyplot as plt 
import matplotlib.tri as tri 
with open('/home/hlias/Desktop/CGSHOP_instances_SA/ortho_20_e2aff192.instance.json', 'r') as f:            # replace for input json file
    input_data = json.load(f)

with open('../output.json', 'r') as f:              # replace for output json file
    output_data = json.load(f)

# original points
points_x = input_data['points_x']
points_y = input_data['points_y']
region_boundary = input_data['region_boundary']
additional_constraints = input_data['additional_constraints']

# steiner points
steiner_points_x = [eval(x) for x in output_data['steiner_points_x']]  # covert to float 
steiner_points_y = [eval(y) for y in output_data['steiner_points_y']]  
edges = output_data['edges']

# merge points
all_points_x = points_x + steiner_points_x
all_points_y = points_y + steiner_points_y

# create figure
fig, ax = plt.subplots(figsize=(10, 10))

# Plot the region boundary
region_points = [(points_x[i], points_y[i]) for i in region_boundary]
region_x, region_y = zip(*region_points)
region_x += (region_x[0],)  
region_y += (region_y[0],)  
ax.plot(region_x, region_y, color='green', label='Region Boundary / Constraints')

# plotting
for constraint in additional_constraints:
    x_coords = [points_x[constraint[0]], points_x[constraint[1]]]
    y_coords = [points_y[constraint[0]], points_y[constraint[1]]]
    ax.plot(x_coords, y_coords, color='green')
for edge in edges:
    x_coords = [all_points_x[edge[0]], all_points_x[edge[1]]]
    y_coords = [all_points_y[edge[0]], all_points_y[edge[1]]]
    ax.plot(x_coords, y_coords, color='blue', alpha=0.6)
ax.scatter(points_x, points_y, color='black', label='Original Points')
ax.scatter(steiner_points_x, steiner_points_y, color='red', label='Steiner Points')

ax.legend()

ax.set_aspect('equal', adjustable='box')
plt.title('Solution CDT')
plt.show()