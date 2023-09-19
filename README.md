# Pixel centroiding using DBSCAN

## Overview

This is a C++ script designed to perform centroiding using the DBSCAN algorithm for clustering given input from a timepix detector. It offers way faster performance compared to Scikit-learn's DBSCAN  method. The script can be called from Python using a function called `run_centroiding`.

Compilation of the C++ script is done automatically.


## Prerequisites

- **C++ Compiler**: You need a C compiler installed on your system to compile and run the C++ script.

## Usage

### From Python

You can call the centroiding script from your Python code using the `run_centroiding` function provided in the `centroider.py` module. Here's an example:

```python
from centroider import run_centroiding

# Example usage with correctiondata (optional)
input_file = "input_data.csv"
output_file = "output_file.csv"
correction_data = "correction_data.csv"  # Optional

run_centroiding(input_file, output_file, correction_data)
```

This function will automatically compile the C++ script if not already done and run it with the given file paths.

### Function Arguments
inputfile (str): Path to the input data CSV file.
outputfile (str): Path to the output centroid CSV file.
correctiondata (str, optional): Path to the correction data CSV file (optional)

# Input File Format
The input data CSV file should have 6 columns with the following headers:

(Index, shot, x, y, tof, tot)

The output centroid CSV file will have the same columns and format as the input file.
 
# Acknowledgments
This project was developed by Laurits Lønskov Sørensen and Simon Fischer-Nielsen and is not affiliated with or endorsed by the owner of Timepix 

Please feel free to report any issues or contribute to this project on GitHub.

