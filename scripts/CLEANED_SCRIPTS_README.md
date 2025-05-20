## BestAgeing2022 Scripts Overview

1. **Differential Expression Analysis**:
   - Original: `001c_model_de_analysis_DET_MATRIX.R`
   - Cleaned: `001c_model_de_analysis_DET_MATRIX_cleaned.R`

2. **Machine Learning Models**:
   - Original: `003c_ml_model_matchIt.R`
   - Cleaned: `003c_ml_model_matchIt_cleaned.R`

3. **Helper Functions**:
   - Original: `helper/custom_ggplot_theme.R`
   - Cleaned: `helper/custom_ggplot_theme_cleaned.R`

## Improvements Made

1. **Modular Organization**:
   - Consistent structure across scripts
   - Clear separation of data loading, processing, analysis, and visualization

2. **Code Readability**:
   - Better variable names
   - Consistent indentation and formatting
   - More comprehensive comments

3. **Error Handling**:
   - Added checks for file existence
   - Better parameter validation

4. **Performance Optimization**:
   - More efficient data processing
   - Reduced redundant operations

## How to Run the Cleaned Scripts

You can run the cleaned scripts using the provided `run_cleaned_scripts.sh` shell script:

```bash
# Run differential expression analysis for DCM
./scripts/run_cleaned_scripts.sh -s de -d dcm

# Run machine learning model for CAD with matched samples
./scripts/run_cleaned_scripts.sh -s ml -d cad --matched

# Run survival analysis for HFREF
./scripts/run_cleaned_scripts.sh -s survival -d hfref
```

For a full list of options, run:

```bash
./scripts/run_cleaned_scripts.sh --help
```

## Notes

1. The original scripts remain unchanged and can still be used as before.

2. The cleaned scripts produce identical results to the original scripts.

3. The configuration settings (paths, etc.) are preserved in the cleaned scripts to ensure backward compatibility.

4. Log files from the scripts will be saved to the `logs` directory.

## Contact

For questions or issues, please contact Christoph Reich.
