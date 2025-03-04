
NebulaPy
========

NebulaPy is a Python library designed to generate line luminosities from PION simulation data. It also facilitates energy binning of stellar atmosphere models, including ATLAS, Potsdam, and CMFGEN, across a wide range of metal abundances. These functionalities are tailored for radiative sources used in PION simulations.

Key Features
------------
- Spectral Energy Distribution (SED) Binning:
  - Supports binned SEDs from models like ATLAS, Potsdam, CMFGEN, and blackbody sources.
- Line Luminosity Calculations:
  - Computes line luminosities for spherically symmetric geometries.
- 2D Cooling Function Maps:
  - Generates 2D maps of cooling functions for astrophysical simulations.
- 2D Line Emissivity Maps:
  - Creates 2D emissivity maps for various line transitions.

Installation and Setup (Linux: Debian/Ubuntu)
---------------------------------------------
1. Clone the repository:
   ```
   git clone https://github.com/arunmathewofficial/NebulaPy.git
   cd NebulaPy
   ```

2. Install dependencies:
   ```
   pip install -r requirements.txt
   ```

3. Download and set up the CHIANTI database:
   - Download the CHIANTI database (~1 GB) from:
     https://download.chiantidatabase.org/CHIANTI_10.1_database.tar.gz
   - Extract it to a directory (~5 GB disk space required):
     ```
     tar -xzf CHIANTI_10.1_database.tar.gz -C ~/MY_CHIANTI_DIRECTORY
     ```
   - Add the following environmental variable to your `.bashrc`:
     ```
     export XUVTOP=$HOME/MY_CHIANTI_DIRECTORY
     ```
     Then, reload your `.bashrc`:
     ```
     source ~/.bashrc
     ```

4. Install the Python-SILO interface:
   - Run the provided script:
     ```
     bash install_silo.sh
     ```
   - Fix the SILO library path in your local distribution:
     - Open the file:
       `${HOME}/.local/venv/lib/python3.11/site-packages/pypion/SiloHeader_data.py`
     - Modify line 18 to append `/lib` to the path.

5. Install NebulaPy:
   ```
   pip install .
   ```

6. To download the database: 
   - Run the following command from the NebulaPy root directory:
     ```bash
     download-database [destination_path]
     ```
   - This requires ~270 MB of additional space.

7. Add environmental variable for NebulaPy Database
   - Add the following environmental variable to your `.bashrc`:
     ```
     export NEBULAPYDB=$HOME/MY_NEBULAPY_DATABASE_DIRECTORY
     ```
     Then, reload your `.bashrc`:
     ```
     source ~/.bashrc
     ```
    
Usage
-----
For usage details, visit the NebulaPy Wiki:
https://github.com/arunmathewofficial/NebulaPy/wiki

Sample Scripts
--------------
Sample scripts demonstrating NebulaPy functionalities can be found in the `NebulaPy/problems` directory.

Support
-------
For bug reports and feature requests, visit the Issues section of the repository:
https://github.com/arunmathewofficial/NebulaPy/issues

Badges
------
License: MIT (https://opensource.org/licenses/MIT)  
PyPI: https://pypi.org/project/NebulaPy/

Author
------
Arun Mathew  
Astronomy & Astrophysics  
Computational and High Energy Astrophysics  
Dublin Institute for Advanced Studies (DIAS), Ireland  

Date: 2024
