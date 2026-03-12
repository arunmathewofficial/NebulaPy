import os
import re
import pandas as pd

input_file = "/Users/tony/Desktop/CWBs-NEMOv1/Post-Processing/Dominant-Lines/wr140_NEMO_d07e13_d2l6n128_DominantLines.txt"
output_dir = "/Users/tony/Desktop/CWBs-NEMOv1/Post-Processing/Dominant-Lines"
output_file = os.path.join(output_dir, "Soft&HardXrayList.txt")
os.makedirs(output_dir, exist_ok=True)

# -------------------------------------------------
# X-ray band limits
# -------------------------------------------------
SOFT_MIN, SOFT_MAX = 0.3, 2.0   # keV
HARD_MIN, HARD_MAX = 2.0, 10.0  # keV

def angstrom_to_keV(wavelength_A):
    return 12.39841984 / wavelength_A

# -------------------------------------------------
# Read file
# -------------------------------------------------
steps = {}
current_step = None
current_time_s = None

with open(input_file, "r") as f:
    for raw_line in f:
        line = raw_line.strip()

        if not line:
            continue

        if line.startswith("Step:"):
            m = re.search(r"Step:\s*(\d+)\s*\|\s*Simulation Time:\s*([0-9.eE+\-]+)", line)
            if m:
                current_step = int(m.group(1))
                current_time_s = float(m.group(2))
                if current_step not in steps:
                    steps[current_step] = {
                        "simulation_time_s": current_time_s,
                        "soft": [],
                        "hard": []
                    }
            continue

        if line.startswith("#"):
            continue

        if "No free-bound emission associated with" in line:
            continue

        parts = line.split()
        if len(parts) < 4:
            continue

        try:
            element = parts[0]
            ion = parts[1]
            wavelength_str = parts[2]
            wavelength_A = float(wavelength_str)
            luminosity = float(parts[3])

            energy_keV = angstrom_to_keV(wavelength_A)
            line_label = f"{element} {ion} {wavelength_str}"

            row = {
                "line": line_label,
                "luminosity": luminosity,
                "energy_keV": energy_keV
            }

            if SOFT_MIN <= energy_keV < SOFT_MAX:
                steps[current_step]["soft"].append(row)
            elif HARD_MIN <= energy_keV <= HARD_MAX:
                steps[current_step]["hard"].append(row)

        except ValueError:
            continue

# -------------------------------------------------
# Sort by luminosity within each band
# -------------------------------------------------
for step in steps:
    steps[step]["soft"].sort(key=lambda x: x["luminosity"], reverse=True)
    steps[step]["hard"].sort(key=lambda x: x["luminosity"], reverse=True)

# -------------------------------------------------
# Write one combined output file
# -------------------------------------------------
with open(output_file, "w") as out:
    for step in sorted(steps.keys()):
        soft_rows = steps[step]["soft"]
        hard_rows = steps[step]["hard"]

        if not soft_rows and not hard_rows:
            continue

        sim_time = steps[step]["simulation_time_s"]

        out.write(f"Step: {step} | Simulation Time: {sim_time:.6e} (s)\n\n")

        if soft_rows:
            out.write("Soft X-ray lines\n")
            out.write(f"{'Line':<24} {'Luminosity':>15} {'Energy (keV)':>15}\n")
            out.write("-" * 60 + "\n")
            for row in soft_rows:
                out.write(
                    f"{row['line']:<24} "
                    f"{row['luminosity']:>15.6e} "
                    f"{row['energy_keV']:>15.6f}\n"
                )
            out.write("\n")

        if hard_rows:
            out.write("Hard X-ray lines\n")
            out.write(f"{'Line':<24} {'Luminosity':>15} {'Energy (keV)':>15}\n")
            out.write("-" * 60 + "\n")
            for row in hard_rows:
                out.write(
                    f"{row['line']:<24} "
                    f"{row['luminosity']:>15.6e} "
                    f"{row['energy_keV']:>15.6f}\n"
                )
            out.write("\n")

        out.write("=" * 72 + "\n\n")

# -------------------------------------------------
# Display on screen
# -------------------------------------------------
for step in sorted(steps.keys()):
    soft_rows = steps[step]["soft"]
    hard_rows = steps[step]["hard"]

    if not soft_rows and not hard_rows:
        continue

    sim_time = steps[step]["simulation_time_s"]

    print(f"Step: {step} | Simulation Time: {sim_time:.6e} (s)")
    print()

    if soft_rows:
        print("Soft X-ray lines")
        print(f"{'Line':<24} {'Luminosity':>15} {'Energy (keV)':>15}")
        print("-" * 60)
        for row in soft_rows:
            print(
                f"{row['line']:<24} "
                f"{row['luminosity']:>15.6e} "
                f"{row['energy_keV']:>15.6f}"
            )
        print()

    if hard_rows:
        print("Hard X-ray lines")
        print(f"{'Line':<24} {'Luminosity':>15} {'Energy (keV)':>15}")
        print("-" * 60)
        for row in hard_rows:
            print(
                f"{row['line']:<24} "
                f"{row['luminosity']:>15.6e} "
                f"{row['energy_keV']:>15.6f}"
            )
        print()

    print("=" * 72)

print(f"\nOutput written to:\n{output_file}")