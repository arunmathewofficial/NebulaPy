import NebulaPy.src as nebula
import pandas
import os

dummy_temperature = [3e+4]
dummy_ne = [1e+3]

output_path = '/home/tony/Desktop/multi-ion-bowshock/NebulaPy/lines-transitions'

nebula_pyneb = nebula.pyneb(pion_ion='H', temperature=dummy_temperature, ne=dummy_ne, verbose=True)

spectroscopic_name = nebula_pyneb.Spectroscopic

lines = nebula_pyneb.get_allLines()
# Convert to a pandas Series for vertical display
series = pandas.Series(lines)
# Save the Series to a text file
filename = os.path.join(output_path, spectroscopic_name.replace(" ", "") + '_recomb_lines.txt')
print(f" saving {spectroscopic_name} recombination lines to file ")
with open(filename, 'w') as f:
    f.write(series.to_string())
