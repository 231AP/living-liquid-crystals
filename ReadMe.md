# LLC-ABP Simulation and Data Analysis Project

This project provides a GPU-accelerated simulation program for active matter, specifically designed for studying the LLC-ABP system. It offers comprehensive capabilities for simulation, data export, visualization, and analysis.

## Features

- **High-Performance Computing**: Leverages CUDA for accelerated simulation on NVIDIA GPUs.
- **Comprehensive Visualization**: Automatically generates images and videos of the simulation process.
- **Flexible Parameterization**: Easily modifiable simulation parameters to suit various research scenarios.
- **Integrated Analysis**: Includes dedicated tools for post-processing and data analysis.

## Requirements

### Hardware
- A computer or cloud server with an **NVIDIA GPU**.

### Software
- The corresponding version of the **NVIDIA CUDA Toolkit** installed.
- **Python 3** (It is recommended to use the latest stable release).

## Installation & Setup

1. **Get the Code**
   Clone this repository to your local machine.

```bash
LLC-ABP
```
2. **Install Python Dependencies**
Install the required Python packages using `pip`.

## Quick Start

1. Navigate to the simulation directory.

2. Run the main script to start the simulation.

3. After the simulation completes, the results will be saved automatically:
- **Data Files**: Located in the `../data/` directory.
- **Images & Videos**: Located in the `../photo_video/` directory.

## Customizing the Simulation

All simulation parameters are defined within the `LLCparams` class in the `LLC-ABP.py` file. You can customize your simulation by modifying the attributes of the `LLCparams` class either directly in `LLC-ABP.py` or via `main.py`.

## Data Analysis

The `analysis` folder in the project root contains helper scripts for data analysis, enabling deeper investigation into the generated simulation data.

## Support & Contact

If you encounter any problems or have suggestions for improvements, please feel free to reach out:
- **Email**: 770395058@qq.com
- **QQ**: 770395058 (Please state your purpose when adding).
