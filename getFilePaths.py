import os
from collections import defaultdict

root_dir = os.path.expanduser("~/Software/soot-dem/dist/soot-dem/example/linux-tutorial/restructuring_simulations")

particles_by_agg_frac = defaultdict(list)

for dirpath, dirnames, filenames in os.walk(root_dir):
    for filename in filenames:
        if filename.startswith("particles_") and filename.endswith(".vtk"):
            parts = dirpath.split(os.sep)
            if len(parts) >= 3 and parts[-3].startswith("aggregate_") and parts[-2].startswith("frac_"):
                agg = parts[-3].split("_")[-1]
                frac = parts[-2].split("_")[-1]
                key = (f"aggregate_{agg}", f"frac_{frac}")
                full_path = os.path.join(dirpath, filename)
                particles_by_agg_frac[key].append(full_path)

for key, files in particles_by_agg_frac.items():
    print(f"{key}: {len(files)} particle files")
    for f in files:
        print("  ", f)
