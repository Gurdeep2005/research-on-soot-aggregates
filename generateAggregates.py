import subprocess
import os
import shutil
import random
import xml.etree.ElementTree as ET
import re

def update_seed_in_config(config_path, new_seed):
    tree = ET.parse(config_path)
    root = tree.getroot()
    for elem in root.findall(".//let"):
        if elem.attrib.get("id") == "rng_seed":
            elem.text = str(new_seed)
            break
    tree.write(config_path)

numberOfAggreagtesToGenerate = 95
bufferAmount = 10
numberOfParticles = 120
randomSeeds = []
path = "/home/gurdeep/Desktop/aggregate_generation/"
aggregationExecutablePath = "/home/gurdeep/Software/soot-dem/build/aggregation"

k0 = 1.3
r = 1.4e-8
Df = 1.8

os.chdir(path)
i = 5
while i < numberOfAggreagtesToGenerate:
    folder_name = f"aggregate_{i}"
    os.makedirs(folder_name)
    shutil.copy("config", folder_name)
    os.chdir(os.path.join(path, folder_name))

    # get a new unused random seed
    randomSeed = random.randint(0, 100000)
    while randomSeed in randomSeeds:
        randomSeed = random.randint(0, 100000)
    randomSeeds.append(randomSeed)

    update_seed_in_config("config", randomSeed)

    result = subprocess.run([aggregationExecutablePath, "config"], capture_output=True, text=True)

    # Search for the radius of gyration line
    match = re.search(r"radius of gyration:\s*([0-9.eE+-]+)", result.stdout)
    if match:
        radius_of_gyration = float(match.group(1))
        os.chdir(path)
        if abs(numberOfParticles - k0*(radius_of_gyration/r)**Df) <= bufferAmount:
            print(f"aggregate_{i} generation succesfull")
            i += 1  # Only move on if simulation was successful and meets equasion
        else:
            shutil.rmtree(folder_name)
            print("Didn't meet equasion requirements! Retrying")
            print("Radius of Gyration: ", radius_of_gyration, "Equasion result: ", abs(numberOfParticles - k0*(radius_of_gyration/r)**Df))
    else:
        print(f"[Aggregate {i}] Failed to extract Radius of Gyration. Retrying...")
        os.chdir(path)
        shutil.rmtree(folder_name)
