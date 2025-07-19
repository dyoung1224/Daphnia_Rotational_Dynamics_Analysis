# Daphnia_Rotational_Dynamics_Analysis

This repository contains a MATLAB script for analyzing the cumulative angular dynamics of tracked Daphnia motion data.

The script processes `.csv` files containing time-series motion data for each individual organism (extracted from TREX `.npz` files), filters them using a velocity threshold (`v*`), and calculates the cumulative angle of rotation over time. It then averages the angular motion across all Daphnia in the experiment and visualizes the result.

## Research Background

This code was developed for an undergraduate research project at **CUNY Queens College** under the supervision of **Dr. Sebastian Alvarado**. The goal of the project is to investigate **behavioral and rotational trends** in Daphnia in response to environmental stimuli, using computational analysis and motion tracking.

This script represents the second step of the analysis pipeline, following the extraction of `.csv` motion data from the raw TREX `.npz` tracking files.

## Features

- Cleans and filters velocity data using a physics-based threshold (`v*`)
- Computes cumulative change in angular direction per Daphnia
- Interpolates angle data to a unified time grid
- Computes and visualizes the average angular behavior across all tracked individuals
- Plots include:
  - Number of Daphnia contributing to analysis over time
  - Average cumulative angle over time

## Input

- A folder of `.csv` files (one per organism), each containing:
  - X, Y position data
  - Time stamps
  - Frame-per-second (fps) information

## Output

- Two plots:
  - Daphnia count over time
  - Average cumulative angular motion over time

## Use Case

This script is designed for analyzing **rotational bias** or **collective angular behavior** in Daphnia motion experiments. It is particularly useful for identifying group behavioral tendencies in response to controlled experimental conditions.

## Dependencies

- MATLAB (tested on R2022b)
- No additional toolboxes required
