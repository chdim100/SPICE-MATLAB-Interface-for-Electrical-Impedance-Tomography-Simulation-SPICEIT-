# SPICE-MATLAB Interface for Electrical Impedance Tomography Simulation (SPICEIT)

Instructions: The main script files are testbench_Tank.m and testbench_Thorax.m

To run them, SPICE transient data must be acquired, using the configurations provided in the Circuit folder (.asc LT SPICE files). 
Make sure that you have all the essential models installed at LT SPICE. 

To create an RLC-equivalent model, use the Set_Thorax.mat function. 

Raw measurement data will be available online soon. 
Please make sure that any data follows the name documentation/pattern defined in the code, so it can be loaded. 

The NETGEN software and EIDORS library are required. 
