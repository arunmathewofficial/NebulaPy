import os
import re
import pandas as pd

input_file = "/Users/tony/Desktop/CWBs-NEMOv1/Post-Processing/Dominant-Lines/wr140_NEMO_d07e13_d2l6n128_DominantLines.txt"
output_dir = "/Users/tony/Desktop/CWBs-NEMOv1/Post-Processing/Dominant-Lines"
output_file = os.path.join(output_dir, "Soft&HardXrayList.txt")
os.makedirs(output_dir, exist_ok=True)

# -------------------------------------------------
# Energy conversion
# -------------------------------------------------
def angstrom_to_keV(w):
    return 12.39841984 / w


SOFT_MIN, SOFT_MAX = 0.3, 2.0
HARD_MIN, HARD_MAX = 2.0, 10.0

# -------------------------------------------------
# Storage
# -------------------------------------------------
steps = {}
current_step = None
current_time = None

# -------------------------------------------------
# Read file
# -------------------------------------------------
with open(input_file) as f:
    for line in f:

        line = line.strip()

        if not line:
            continue

        if line.startswith("Step:"):

            m = re.search(r"Step:\s*(\d+)\s*\|\s*Simulation Time:\s*([0-9.eE+\-]+)", line)

            if m:
                current_step = int(m.group(1))
                current_time = float(m.group(2))

                steps[current_step] = {
                    "time": current_time,
                    "soft": [],
                    "hard": [],
                    "other": []
                }

            continue

        if line.startswith("#"):
            continue

        parts = line.split()

        if len(parts) < 4:
            continue

        try:
            element = parts[0]
            ion = parts[1]
            wavelength_str = parts[2]
            luminosity = float(parts[3])

            wavelength = float(wavelength_str)

            energy = angstrom_to_keV(wavelength)

            line_name = f"{element} {ion} {wavelength_str}"

            row = (line_name, luminosity, energy)

            if SOFT_MIN <= energy < SOFT_MAX:
                steps[current_step]["soft"].append(row)

            elif HARD_MIN <= energy <= HARD_MAX:
                steps[current_step]["hard"].append(row)

            else:
                steps[current_step]["other"].append(row)

        except:
            continue

# -------------------------------------------------
# Sort lines by luminosity
# -------------------------------------------------
for step in steps:

    for band in ["soft", "hard", "other"]:
        steps[step][band].sort(key=lambda x: x[1], reverse=True)

# -------------------------------------------------
# Write ONE output file
# -------------------------------------------------
with open(output_file, "w") as out:
    for step in sorted(steps):

        time = steps[step]["time"]

        header = f"Step: {step} | Simulation Time: {time:.6e} (s)"

        print(header)
        out.write(header + "\n\n")


        def print_block(title, data):

            if not data:
                return

            print(title)
            out.write(title + "\n")

            col = f"{'Line':<24}{'Luminosity':>15}{'Energy (keV)':>15}"

            print(col)
            out.write(col + "\n")

            line = "-" * 55

            print(line)
            out.write(line + "\n")

            for line_name, lum, energy in data:
                row = f"{line_name:<24}{lum:>15.6e}{energy:>15.4f}"

                print(row)
                out.write(row + "\n")

            print()
            out.write("\n")


        print_block("Soft X-ray lines", steps[step]["soft"])
        print_block("Hard X-ray lines", steps[step]["hard"])
        print_block("Other bright lines", steps[step]["other"])

        sep = "=" * 70

        print(sep)
        print()

        out.write(sep + "\n\n")

print("Output file:")
print(output_file)
