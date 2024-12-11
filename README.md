# Crocodile
Cross-shore Coastal Diffusion Long-term Evolution model

# Overview
Crocodile is a diffusion-based cross-shore numerical model designed to simulate the decadal-scale morphological evolution of nourished sandy coasts. It is a simple, robust and computationally efficient tool that provides insights into key coastal indicators like profile volume, coastline position, and beach width under various nourishment strategies.
The model is particularly suited for analysing:

- Long-term depth-dependent cross-shore sand redistribution in sandy coastal profiles responding to nourishment and sea level rise.

- Effects of beach, shoreface and mega nourishment strategies on coastal dynamics. It enables simulating both reactive (e.g. hold-the-line) and proactive (predefined nourishment timing and volume) nourishment strategies.

- Insights into profile steepening, nourishment lifetimes, and volume requirements.

This model has been applied to case studies along the Dutch coastline, demonstrating its utility in predicting long-term morphological trends and optimizing nourishment strategies. 

# How to Use
Follow the steps below to set up and use the model:

## Install Dependencies
   
Follow the steps below to create a conda environment with the dependencies:

1. Create a conda environment using the provided environment.yml file. This will download and install the dependencies required to run the model.

```
conda env create -f environemnt.yml
```

2. Activate the conda environment created above.
```
conda activate crocodile
```



## Preparation to run the model
   
1. Update the Project Path

    - Open the Crocodile_input.py file.
    
    - Modify the project_directory variable to reflect the path to your project folder. Example: project_directory = "/path/to/your/project"

2. Create the Required Directory

    - In the root of your project folder, create an empty directory named `0`:

```
mkdir 0
```

3. Grant write permissions to the directory:

```
chmod u+w 0
```

## Run the Model

To execute the model, follow these steps:

- Navigate to the project folder in your terminal.

- Run the main script:

```
python Crocodile_input.py
```

The output will include updated coastal profile simulations and visualizations of morphological changes.

# File Structure

Key Files

- Crocodile_input.py: Configuration and setup for the project.

- Crocodile_model.py: Core script for running the simulation.

- Analysis_functions.py: Analytical utilities for processing simulation outputs.

- Case_study_functions.py: Tools for handling case-specific data.

- Profile_functions.py: Functions for generating and manipulating coastal profiles.

- Plot_CS_profiles.py: Visualization tools for cross-shore profile analysis.

- plot_scenario_snaps.py: Scenario-specific plotting.

- transects.py: Utilities for managing JARKUS transect data.


# Input Data
Input data such as bathymetry, nourishment configurations, and transect details should be organized in subfolders specified in Crocodile_input.py.

# Output Data
Simulation results, including netCDF files and visual plots, will be saved in the respective scenario folders.

# Acknowledgements 
This model was developed as part of the research project CScape. For further details, refer to Kettler et al. (2024): https://doi.org/10.1016/j.coastaleng.2024.104491.

