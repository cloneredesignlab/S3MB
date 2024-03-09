Overview
Manager.m orchestrates the GBM growth simulation and parameter optimization process, predicting future tumor growth and treatment responses. It initializes simulations using AllParams_transition2real.xlsx as the starting parameter set.

Key Components

Data Files:

AllParams_transition2real.xlsx: Used initially to start the simulation with a baseline parameter set. Represents patient-independent parameters.

AllParams.xlsx: Generated or utilized after fitting S3MB simulations to clinical data by minimizing the cost function, incorporating optimized parameters for specific patient data.

Special Flags:

GROWELSEWHERE: A boolean flag within Manager.m that, when set to true, flips the tumor growth and resection cavity data matrices. This is used to simulate or adjust the orientation of tumor growth potentially in different brain areas or to accommodate specific anatomical considerations in the simulations.

Simulation Process

Initialization with AllParams_transition2real.xlsx: The script begins with this parameter file to configure the initial conditions for the simulation.
Parameter Optimization: Through simulations and using cost.m, it adjusts parameters to align simulated tumor growth with observed clinical data, updating or creating AllParams.xlsx with these optimized parameters.
Prediction and Visualization: Utilizing the optimized parameters, it predicts future tumor growth and responses to treatments, with options to visualize these predictions and compare them against clinical observations.

Running the Simulation
Ensure all prerequisites, including MATLAB and necessary toolboxes, are met. Patient-specific data should be in place, and both initial and optimized parameter files should be correctly referenced within the script. Run Manager.m. 
