import NebulaPy.src as nebula
from NebulaPy.tools import util
from pypion.ReadData import ReadData
import astropy.units as unit

# Set up paths and filenames
silo_dir = '/home/tony/Desktop/NebulaPy/tests/wind-wind-jm'  # Directory containing silo files
filebase = 'e7_WRwind_d1l5n256_v0750'  # Base name of the silo files
output_file = '/home/tony/Desktop/NebulaPy/tests/line_luminosity_OIV25.txt'  # Output file for results

# Set up the ion and line emission parameters
ion_name = 'O3+'  # The ion of interest (Oxygen IV)
mass_O = 2.6567628e-23  # Mass of an Oxygen atom (in grams)
line = 2.589332e+05  # Emission line of interest
print(rf" calculating line luminosity of {ion_name} {line} \AA")


# Batch the silo files according to the time instant
batched_silos = util.batch_silos(silo_dir, filebase)

# Initialize the Pion class from NebulaPy, which handles the simulation data
nebula_pion = nebula.pion(batched_silos, verbose=True)

# Extract all chemistry information from the silo files into a chemistry container
# This uses the first time instant's silo file to initialize
nebula_pion.get_chemistry()

# Initialize spherical grid parameters (e.g., radius, shell volumes)
# This sets up the grid using the first silo file in the batch
nebula_pion.spherical_grid(batched_silos[0])

# Retrieve the radius and shell volumes from the geometry container
radius = nebula_pion.geometry_container['radius']
shell_volume = nebula_pion.geometry_container['shell_volumes']

line_emission = nebula.emissionline(ion_name, verbose=True)  # Initialize the emission line calculation

# Open the output file and write the header
with open(output_file, "w") as file:
    file.write(util.nebula_version()+'\n')
    file.write(f"time (kyr)\tLine {line} Luminosity (erg/s)\n")

# Loop over each time instant in the batched silo files
for silo_instant in batched_silos:
    # Read the data from the current silo file
    dataio = ReadData(silo_instant)
    basic = dataio.get_1Darray('Density')  # Retrieve basic simulation data, such as density
    time = (basic['sim_time'] * unit.s).to(unit.kyr)  # Convert simulation time to kiloyears
    dataio.close()  # Close the data file

    # Print the current time instant
    print(f" ---------------------------")
    print(f" time {time:.6e}")

    # Extract necessary physical parameters for the current time instant
    density = nebula_pion.get_parameter('Density', silo_instant)  # Retrieve density
    temperature = nebula_pion.get_parameter('Temperature', silo_instant)  # Retrieve temperature
    ion_massfrac = nebula_pion.get_ion(ion_name, silo_instant)  # Retrieve ion mass fraction

    # Calculate ion number density (number of ions per unit volume)
    ion_num_density = ion_massfrac * density / mass_O

    # Retrieve the electron number density
    ne = nebula_pion.get_ne(silo_instant)

    # Calculate the line luminosity for the specific emission line
    line_emission.lineluminosity_spherical(
        line=line,
        temperature=temperature,
        ne=ne,
        ns=ion_num_density,
        dV=shell_volume
    )

    # Print the calculated luminosity for the current time instant
    print(f" L_{line} = {line_emission.LineLuminosity['luminosity']} erg/s")

    # Append the calculated luminosity to the output file
    with open(output_file, "a") as file:
        file.write(f"{time.value:.6e}\t{line_emission.LineLuminosity['luminosity']:.6e}\n")

