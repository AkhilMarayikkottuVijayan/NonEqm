# NonEqm
NonEQM UDF Ansys


A sample UDF for defining the vibrational energy state of N2 species as a function of temperature. Simple Harmonic Oscillator moel is used in this case to determine the vib energy state. 

Things to do:

1. Make this routine genertic for all species by making the molecular weight MW and characteristic vib temperature for teh species as an input variable
2. Develop a vib energy source term SV to couple the vib mode with the total energy eqn
3. Use UDS to model the transport of vib energy state in the system
