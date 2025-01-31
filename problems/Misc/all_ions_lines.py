import NebulaPy.src as nebula
import pandas
import os

temperature = [2e+6]
ne = [1e+9]

output_path = '/home/tony/Desktop/Equi_NonEqui/'

nebula_chianti = nebula.chianti(pion_ion='He1+', temperature=temperature, ne=ne, verbose=True)

spectroscopic_name = nebula_chianti.chianti_ion.Spectroscopic

lines = nebula_chianti.get_allLines()
# Convert to a pandas Series for vertical display
series = pandas.Series(lines)
# Save the Series to a text file
filename = os.path.join(output_path, spectroscopic_name + '_lines.txt')
print(f" saving {spectroscopic_name} lines to file ")
with open(filename, 'w') as f:
    f.write(series.to_string())
