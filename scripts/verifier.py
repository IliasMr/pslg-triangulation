"""
=======================================================

Python script from eclass to verify the solution
Input/Output should be changed for every instance

=======================================================
"""

from cgshop2025_pyutils.data_schemas.instance import Cgshop2025Instance
from cgshop2025_pyutils.data_schemas.solution import Cgshop2025Solution
from cgshop2025_pyutils import verify
import json

# input_filepath = "/home/hlias/Desktop/CGSHOP_instances_random/ortho_20_e2aff192.instance.json"
input_filepath = "/home/hlias/Desktop/CGSHOP_instances_SA_L_300/ortho_100_bd1e4a14.instance.json"           # replace for input
output_filepath = "/home/hlias/Desktop/CG:SHOP-solutions-SA_L_300/ortho_100_bd1e4a14.instance_output.json"            # replace for output

with open(input_filepath, 'r') as file:
  data_in = json.load(file)

with open(output_filepath, 'r') as file:
  data_out = json.load(file)


instance = Cgshop2025Instance(
            instance_uid=data_in["instance_uid"],
            num_points=data_in["num_points"],
            points_x=data_in["points_x"],
            points_y=data_in["points_y"],
            region_boundary=data_in["region_boundary"],
            num_constraints=data_in["num_constraints"],
            additional_constraints=data_in["additional_constraints"],
          )


solution = Cgshop2025Solution (
              content_type="CG_SHOP_2025_Solution",
              instance_uid=data_out["instance_uid"],
              steiner_points_x=data_out["steiner_points_x"],
              steiner_points_y=data_out["steiner_points_y"],
              edges=data_out["edges"],
            )

result = verify(instance, solution, strict=True)

if result.num_obtuse_triangles != -1:
  print(f"No. obtuse triangles: {result.num_obtuse_triangles}\nNo. Steiner points: {result.num_steiner_points}")
else:
  print("Errors:")
  for err in result.errors: print(f"- {err}")  