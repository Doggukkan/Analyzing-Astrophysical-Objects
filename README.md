# Black Hole X-ray Emission Analysis Project

## Overview

This project aims to study black holes and their X-ray emissions by analyzing astronomical data through a series of computational steps. Utilizing academic resources, the project will focus on understanding black holes, X-ray emissions, and the capabilities of X-ray observatories. The main goal is to process and visualize the data to identify specific characteristics of the X-ray emissions, particularly in the context of dust-scattering halos around black holes.

## Methodology

### Research Phase

- Conduct thorough literature research on:
  - Black holes
  - X-ray emissions from black holes
  - X-ray observatories

### Data Preparation

- Use the Astropy library in Python to format the FITS file data into a filtered array with [x coordinate, y coordinate, photon energy].
- This array will help in pinpointing specific pixels within the desired image area.

### Data Analysis

- Develop an algorithm to detect pixels within the desired area that match our criteria, focusing on the relevant data points for further examination.

### Visualization

- Transform the filtered data into a histogram to analyze the energy levels within the specified radials using the matplotlib library in Python.
- Utilize the data to plot the radial density profile of a dust-scattering halo, providing insights into the distribution and intensity of X-ray emissions around black holes.

## Tools and Libraries

- **Python:** The primary programming language used for the project.
- **Astropy:** A Python library for astronomy-related data processing and analysis.
- **Matplotlib:** A Python library used for data visualization, including histograms and plots.

## Objective

The project seeks to enhance our understanding of black holes and their surrounding environments by analyzing X-ray emissions. Through sophisticated data analysis and visualization techniques, we aim to uncover patterns and characteristics of dust-scattering halos, contributing to the broader field of astrophysics research.
