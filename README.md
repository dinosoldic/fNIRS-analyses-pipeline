# fNIRS Analyses Pipeline

This repository provides a basic pipeline for fNIRS data analysis using MATLAB. It includes tools to organize data, perform analysis, and visualize results.
Currently, mass-univariate dependent and independent t-tests of NIRS task-related time domain are included.

**Note:** All files should be kept in the same directory since functions are co-dependent.

## Directory Structure

The directory MUST be organized for loadData() as follows:
```
 Study
 ├── Group1
 │   ├── Condition1
 │   │   ├── data1
 │   │   └── dataN
 │   └── Condition2
 │       ├── data1
 │       └── dataN
 ├── Group2
 │   ├── Condition1
 │   │   ├── data1
 │   │   └── dataN
 │   └── Condition2
 │       ├── data1
 │       └── dataN
 └── Group3
     ├── Condition1
     │   ├── data1
     │   └── dataN
     └── Condition2
         ├── data1
         └── dataN
```
## Pipeline Overview

1. **Extract Files (Optional)**
   - Use `extractFiles(parent_root, output_root)` to organize HRF files into subject/group/condition folders.
   - **Note:** This function requires user-specific adjustments and is not general-purpose.

2. **NIRS Analysis**
   - Use `NIRSAnalysis` to perform analyses.
   - This function calls `exportNIRS` internally to save results as `.csv`.
  
3. **Visualization**
   - Use `plotNIRS(results, data)` to visualize fNIRS signals and statistical results.
   - Supports plots with significance markers.

## Example Usage

```matlab
% Optional: organize files
extractFiles("raw_data_folder", "organized_data_folder");

% Run analysis (loadData and exportNIRS called internally)
NIRSAnalysis(data);

% Plot results
plotNIRS(results, data);

% Optional: export results
exportNIRS(results, "export_folder");
