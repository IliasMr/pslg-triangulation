"""
=======================================================

Python script for running multiple instances
Input/Output should be changed manually for every instance like the verifier

=======================================================
"""


import os
import subprocess
import re
from concurrent.futures import ThreadPoolExecutor, as_completed

def process_instance(instance_path, output_folder):
    """Process a single JSON instance."""
    instance_name = os.path.basename(instance_path)
    unique_output_path = os.path.join(output_folder, f"{os.path.splitext(instance_name)[0]}_output.json")

    # run project for every instance 
    process = subprocess.run([
        "/home/hlias/Desktop/coding/7o eksamhno/project/hw3/project3-di/build/./opt_triangulation", "-i", instance_path, "-o", unique_output_path   # replace for run command
    ], capture_output=True, text=True)

    output = process.stdout
    results = []

    # parse the program output and get info 
    obtuse_after_match = re.search(r"Obtuse triangles in cdt after:\s*(\d+)", output)
    energy_match = re.search(r"Energy:\s*([\d.]+)", output)

    if obtuse_after_match and int(obtuse_after_match.group(1)) == 0:
        results.append(f"{instance_name} - Converge rate instead of Energy.\n")
    elif energy_match:
        energy = float(energy_match.group(1))
        results.append(f"{instance_name} - Energy: {energy:.2f}\n")
        return instance_name, output, energy
    else:
        results.append(f"{instance_name} - No Energy or Conv rate found.\n")
        return instance_name, output, None

    return instance_name, output, None

def run_project(input_folder):
    output_file = "/home/hlias/Desktop/SA-results-L_3000-b_0.5.txt"            # replace for output.txt
    output_folder = "/home/hlias/Desktop/instances-solutions-SA_L_3000-b_0.5"                               # replace for outputs folder
    os.makedirs(output_folder, exist_ok=True)  

    total_energy = 0
    energy_count = 0

    files = [f for f in os.listdir(input_folder) if f.endswith(".json")]
    total_files = len(files)

    # using threads for faster execution
    with ThreadPoolExecutor() as executor:
        
        futures = {
            executor.submit(process_instance, os.path.join(input_folder, filename), output_folder): filename
            for filename in files
        }

        with open(output_file, "w") as results:
            for idx, future in enumerate(as_completed(futures), start=1):
                try:
                    instance_name, output, energy = future.result()

                    # append results to the output file
                    results.write("========================\n")
                    results.write(f"{instance_name}\n")
                    results.write(output)
                    results.write("========================\n")

                    if energy is not None:
                        total_energy += energy
                        energy_count += 1

                except Exception as e:
                    results.write(f"Error processing {futures[future]}: {e}\n")

                # progress
                print(f"{idx}/{total_files}")

            # calculate and append the total medium Energy
            if energy_count > 0:
                medium_energy = total_energy / energy_count
                results.write(f"\nTotal Medium Energy: {medium_energy:.2f}\n")
            else:
                results.write("\nNo Energy values were found to calculate Medium Energy.\n")

if __name__ == "__main__":
    input_folder = "/home/hlias/Desktop/CGSHOP_instances_SA_L_300_b_0.5"   # replace for input folder
    run_project(input_folder)
