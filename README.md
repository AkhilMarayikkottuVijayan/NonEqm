# NonEqm
NonEQM UDF Ansys


A sample UDF for defining the vibrational energy state of N2 species as a function of temperature. Simple Harmonic Oscillator moel is used in this case to determine the vib energy state. 

Things to do:

1. Make this routine genertic for all species by making the molecular weight MW and characteristic vib temperature for the species as an input variable

2. Develop a species flux function using UDFs

3. Equilibrium constant formula to be added to "SourceReac.c" code to compute Keq

4. Validation case : (Hornung et al. 1972 )
