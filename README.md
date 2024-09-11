# BioRT-HBV

BioRT-HBV is a reactive transport module that works with the semi-distributed rainfall-runoff model [HBV Light](https://www.geo.uzh.ch/en/units/h2k/Services/HBV-Model.html).
BioRT-HBV takes HBV Light simulated hydrologic states and fluxes to drive reactive transport processes.
The reaction processes are the same as in [BioRT-Flux-PIHM](https://github.com/PSUmodeling/MM-PIHM).

BioRT-HBV is open source software licensed under the MIT License.
All bug reports and feature requests should be submitted using the Issues page.

## Usage

The following guide applies to UNIX (include MacOS) systems.

### Installing BioRT-HBV

#### Linux/MacOS
After downloading the BioRT-HBV source code, go into the `BioRT-HBV` directory.
Use the command

```shell
make all
```

to install CVODE library (required by BioRT-HBV) and compile BioRT-HBV executable.

If you already have CVODE v2.9.0 installed, you can edit the Makefile and point `CVODE_PATH` to your `CVODE` directory, and use

```shell
make biort
```

to compile BioRT-HBV.

When installation succeeds, you should see a `biort` executable in your `BioRT-HBV` directory.

#### Windows
To build BioRT-HBV for Windows, ensure that you have a C compiler and CMake installed on your local machine. Then in a Powershell terminal, run:
```shell
mkdir build
cd build
cmake .. 
cmake --build . --config Release
```

If successfully compiled, the executable `biort.exe` should be in the parent directory.

###windows - without compiling
To directly run the BioRT-HBV model on your system, unzip "Biort_HBV_v1_windows_exe.zip". There will be "biort.exe" file and "input" folder inside "Biort_HBV_v1_windows_exe" folder. 
Navigate to the folder "Biort_HBV_v1_windows_exe" on your command prompt and you can run the model using the command "biort xxx", where xxx is a subfolder of "input" folder containing the input files for the study site user wants to model.

### Preparing input files

The BioRT-HBV model requires six input files, including HBV Light model input and output files.

#### Input files from HBV Light

Before running the BioRT-HBV model, users should first run the HBV Light model for their watersheds of interest.
The parameters used for the HBV Light simulation should be saved.
After HBV Light simulations, create a sub-directory in the `BioRT-HBV` `input` directory.
White spaces should be avoided when naming the sub-directory.
Then put the saved parameter file and `Results.txt` to the created sub-directory.
The parameter file should be renamed to `Parameter.xml`.

Four additional files, `cdbs.txt`, `chem.txt`, `cini.txt`, and `soil.txt` are required for the BioRT-HBV model.

#### Chemical database file

The chemical database, `cdbs.txt`, has the exactly same format as CrunchFlow and BioRT-Flux-PIHM database files.
For a complete description on the CrunchFlow's database file, please refer to CrunchFlow user's manual.

#### Chemical model control file

The chemical model control file, `chem.txt`, follows the reactive transport code CrunchFlow input file structure.

The first block provides parameters to control the reactive transport module:

`RECYCLE`:
Number of times to recycle HBV Light results as forcing for spin-up.
The model will run `RECYCLE + 1` times.

`ACTIVITY`:
Set to 0 to disable activity model.
Set to 1 to use the Debye-Huckel equation.

`RELMIN`:
Mineral volume fraction mode.
Set to 0 if you input absolute volume fraction.
Set to 1 if you input relative fraction (to total mineral volume fraction).

`TRANSPORT_ONLY`:
Transport only mode switch.
Set to 0 to turn on reactions.
Set to 1 to run transport only mode with reactions skipped.

`CEMENTATION`:
Cementation factor.

`TEMPERATURE`:
Field temperature in degree celsius.

`SW_THRESHOLD` and `SW_EXP`:
Threshold and exponent parameters in soil moisture control function.
Set `SW_THRESHOLD` to 1 to use increase behavior.

`Q10`:
*Q*<sub>10</sub> parameter for reactions.

The `PRIMARY_SPECIES` block lists all primary species to be simulated.
The `SECONDARY_SPECIES` block lists all secondary species the users would like to track.
The `MINERAL_KINETICS` block lists all kinetic mineral reactions.
This block follows a similar format as the corresponding block in CrunchFlow input file.

#### Chemical initial condition file

The chemical initial condition file (`cini.txt`) provides precipitation concentrations (`PRECIPITATION`) and the initial concentrations in the upper (`UZ`) and lower (`LZ`) zones.
The total concentrations of all primary species should be provided in the `cini.txt` file.
The default unit for aqueous concentrations is mol&nbsp;L<sup>−1</sup>,
for volume fractions is m<sup>3</sup>&nbsp;m<sup>−3</sup>,
for specific surface area is m<sup>2</sup>&nbsp;g<sup>−1</sup>,
for surface site density is eq&nbsp;g<sup>−1</sup>.

#### Soil parameter file

The soil parameter file (`soil.txt`) adds necessary soil parameters to BioRT simulations.

`POROSITY_UZ` and `POROSITY_LZ`:
Porosity of upper and lower zones (m<sup>3</sup>&nbsp;m<sup>-3</sup>);

`RES_UZ` and `RES_LZ`:
Residual water storages in upper and lower zones (mm).
This storage is the lowest allowed storages in BioRT-HBV simulations.
The residual storages are added to simulated storages at every time step.

`D_UZ` and `D_LZ`:
Depths of upper and lower zones (mm).
The depths are used for the model to calculate degrees of saturation.

### Running BioRT-HBV

Now you can run BioRT-HBV models using:

```shell
$ ./biort [-vV] <input directory>
```

where `<input directory>` is the input directory name of your choice,
and `[-v]` is an optional parameter to turn on verbose mode.

The `-V` parameter will display model version.
Note that model will quit after displaying the version information.
No simulation will be performed when using the `-V` parameter.

### Output file

After running the simulation, an output file containing the concentrations of all species and reaction rates of minerals at each time step is generated in your `BioRT-HBV/output` directory named `<input directory>_results_<simulation_time>.txt`.

The header line shows the names of all species, followed by a suffix indicating locations (and rates).

`_inf`:
Infiltration concentrations (mol L<sup>-1</sup>).

`_UZ`:
Upper zone concentrations (mol L<sup>-1</sup>).

`_LZ`:
Lower zone concentrations (mol L<sup>-1</sup>).

`_riv`:
Stream concentrations (mol L<sup>-1</sup>).

`_rate_UZ`:
Upper zone mineral reaction rates (mol m<sup>-2</sup> day<sup>-1</sup>)

`_rate_LZ`:
Lower zone mineral reaction rates (mol m<sup>-2</sup> day<sup>-1</sup>)
