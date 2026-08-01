# Nanopore Signal Analysis and Z-score-based Detection

[![DOI](https://zenodo.org/badge/1319110691.svg)](https://doi.org/10.5281/zenodo.21737767)

## Overview

This repository contains the custom downstream analysis scripts used for nanopore signal comparison, Z-score-based detection, sequencing-error analysis, statistical testing, classification-performance evaluation, and figure generation.

The repository provides scripts for:

* comparison and visualization of reference-aligned nanopore signal levels;
* position-specific statistical analysis of signal differences;
* Z-score-based classification of control and dA-AL-II reads;
* evaluation of multiple Z-score thresholds;
* confusion-matrix and receiver operating characteristic (ROC) analysis;
* simulated mixture-proportion analysis;
* sequencing-error, odds-ratio, and Fisher's exact-test analysis;
* generation of figures and tabular results used in the study.

The signal-analysis scripts use prepared, normalized, reference-aligned signal matrices as input. The corresponding raw signal files, sequencing reads, reference sequence, prepared signal matrices, and processed analysis results are publicly available in the associated Zenodo dataset:

https://doi.org/10.5281/zenodo.21740475

## Repository contents

### `plot_signal_compare.py`

Compares reference-aligned signal traces between control and dA-AL-II samples for QT and ONT datasets.

The script generates:

* whole-region signal-comparison plots;
* local signal-comparison plots;
* position-specific kernel-density plots.

A median filter with a kernel size of 5 is applied when visualizing aligned signal traces.

### `zscore_read_classifier.py`

Performs position-specific Z-score-based detection.

The script:

* uses the first half of the control reads to calculate the background mean and standard deviation at each reference position;
* uses the second half of the control reads as a held-out control evaluation set;
* calculates the absolute Z-score of each control and dA-AL-II read;
* evaluates Z-score thresholds of 1, 2, and 3;
* calculates false-positive rate, true-positive rate, precision, and F1 score;
* generates normalized confusion matrices;
* generates ROC curves and calculates the area under the curve;
* evaluates predicted positive proportions using simulated mixtures of control and dA-AL-II reads.

The current implementation uses a Z-score threshold of 2.0 for the mixture-proportion analysis.

### `dis_plot.py`

Performs site-level statistical comparison of signal levels between control and dA-AL-II samples.

The script:

* calculates mean signal levels at each analyzed position;
* calculates fold changes;
* performs Welch's two-sample t-tests;
* calculates standardized signal deviations relative to the control distribution;
* generates a radar plot highlighting the target position.

The default platform selected in the current script is QT.

### `generate_odds.py`

Performs position-specific sequencing-error analysis using aligned BAM files and a reference FASTA file.

The script:

* calculates position-specific mismatch and deletion counts;
* calculates combined error rates based on the currently implemented mismatch and deletion counts;
* summarizes base-quality scores;
* compares dA-AL-II and control samples;
* calculates odds ratios using a pseudocount of 0.5;
* performs one-sided Fisher's exact tests;
* generates odds-ratio, error-rate radar, and quality-comparison plots;
* exports position-specific results as a CSV file.

The current script processes a maximum of 10,000 reads per BAM file by default and retains positions with at least 30 reads in both groups for the main comparison.

## Repository structure

```text
nanopore-zscore-analysis/
├── README.md
├── UPLOAD_GUIDE_CN.md
├── LICENSE
├── requirements.txt
├── environment.yml
├── .gitignore
├── dis_plot.py
├── generate_odds.py
├── plot_signal_compare.py
├── zscore_read_classifier.py
├── data/                      # Create locally for input files
└── fig/                       # Created automatically for output figures
```

The `data/` directory is not included in the GitHub repository because it contains large research-data files.

The complete data archive is publicly available in Zenodo:

https://doi.org/10.5281/zenodo.21740475

The `fig/` subdirectories are created automatically when the scripts are run.

## Reproducibility scope and input data

This repository provides the custom downstream analyses used for nanopore signal comparison, statistical testing, Z-score-based detection, classification-performance evaluation, sequencing-error analysis, and figure generation.

The signal-analysis scripts use prepared, normalized, reference-aligned signal matrices as input. The sequencing-error analysis uses basecalled and reference-aligned BAM files together with a custom reference sequence.

The required data files are publicly available in the associated Zenodo dataset:

https://doi.org/10.5281/zenodo.21740475

The Zenodo dataset contains:

* raw Oxford Nanopore FAST5 signal files;
* raw Qitan Technology H5 signal files;
* basecalled and reference-aligned BAM files;
* BAM index files;
* the custom reference sequence;
* prepared reference-aligned signal matrices;
* processed sequencing-error and odds-ratio results;
* a data README describing the deposited files.

This repository focuses on the downstream analytical workflow and does not provide an end-to-end implementation of platform-specific upstream signal processing.

## Software requirements

A compatible Python 3 environment is required.

The main Python dependencies include:

* NumPy;
* pandas;
* SciPy;
* Matplotlib;
* seaborn;
* scikit-learn;
* pysam.

Install the required packages using:

```bash
pip install -r requirements.txt
```

Alternatively, create the provided Conda environment:

```bash
conda env create -f environment.yml
conda activate nanopore-zscore-analysis
```

Using the provided environment specification is recommended.

## Obtaining the input data

Download the associated dataset from Zenodo:

https://doi.org/10.5281/zenodo.21740475

Create a directory named `data` in the repository root:

```bash
mkdir data
```

Place the required files in the `data/` directory before running the analyses.

The expected local structure is:

```text
nanopore-zscore-analysis/
├── data/
│   ├── QT_control_align_res.pkl
│   ├── QT_dAAL_align_res.pkl
│   ├── ONT_control_align_res.pkl
│   ├── ONT_dAAL_align_res.pkl
│   ├── ONT_control.bam
│   ├── ONT_control.bam.bai
│   ├── ONT_dAAL.bam
│   ├── ONT_dAAL.bam.bai
│   ├── QT_control.bam
│   ├── QT_control.bam.bai
│   ├── QT_dAAL.bam
│   ├── QT_dAAL.bam.bai
│   └── ref.fasta
└── analysis scripts
```

Only the files required for a particular analysis need to be placed in the directory.

## Expected signal-matrix inputs

The signal-comparison, DIS, and Z-score scripts use the following prepared signal-matrix files:

```text
data/QT_control_align_res.pkl
data/QT_dAAL_align_res.pkl
data/ONT_control_align_res.pkl
data/ONT_dAAL_align_res.pkl
```

Each pickle file must contain a Python dictionary with a key named:

```python
"sigs"
```

The value associated with `sigs` must be a two-dimensional signal matrix in which:

* each row represents an individual nanopore read;
* each column represents a reference-aligned signal measurement;
* six signal measurements are used per reference position.

The compatible signal matrices must have already undergone the required upstream signal processing, normalization, reference alignment, and resampling.

The prepared signal matrices used by these scripts are publicly available in the associated Zenodo dataset:

https://doi.org/10.5281/zenodo.21740475

## Expected sequencing-analysis inputs

The sequencing-error analysis uses:

```text
data/ONT_control.bam
data/ONT_dAAL.bam
data/ref.fasta
```

Indexed BAM files are required because the script accesses aligned reads by reference sequence:

```text
data/ONT_control.bam.bai
data/ONT_dAAL.bam.bai
```

The BAM files must be aligned to the reference sequence contained in `ref.fasta`.

The current script analyzes the first reference sequence in `ref.fasta`. The reference name used in the BAM files must therefore match the first reference name in the FASTA file.

The corresponding BAM files, BAM index files, and reference sequence are available in the associated Zenodo dataset.

## Running the analyses

Run the scripts from the root directory of the repository.

### Signal comparison and visualization

```bash
python plot_signal_compare.py
```

This script processes both QT and ONT aligned signal matrices and saves results under:

```text
fig/QT/
fig/ONT/
```

### Z-score-based detection

```bash
python zscore_read_classifier.py
```

This script processes both QT and ONT aligned signal matrices and generates:

* Z-score threshold plots;
* normalized confusion matrices;
* ROC curves;
* simulated mixture-proportion heatmaps.

Results are saved under:

```text
fig/QT/
fig/ONT/
```

### Site-level signal-difference analysis

```bash
python dis_plot.py
```

The platform is selected using the `platform` variable inside the script.

The generated radar plot is saved under the corresponding platform directory:

```text
fig/QT/
```

or:

```text
fig/ONT/
```

### Sequencing-error and odds-ratio analysis

```bash
python generate_odds.py
```

The platform is selected using the `platform` variable inside the script.

The script generates:

* a position-specific odds-ratio table;
* an odds-ratio plot;
* an error-rate radar plot;
* a quality-score comparison plot.

The CSV result is saved under:

```text
data/
```

The figures are saved under the corresponding platform directory:

```text
fig/QT/
```

or:

```text
fig/ONT/
```

## Analysis configuration

Several analysis settings are currently defined directly inside the Python scripts.

These settings include:

* sequencing platform (`QT` or `ONT`);
* input filenames;
* analyzed reference-position range;
* reference-coordinate offset;
* signal-resampling factor;
* Z-score thresholds;
* Z-score threshold used for mixture analysis;
* minimum read depth;
* maximum number of BAM reads;
* target position;
* statistical significance threshold;
* output directories.

The current signal-analysis scripts use six aligned signal measurements per reference position.

The principal highlighted target position is reference position 506.

Before applying the scripts to another dataset, review the configuration variables inside the relevant script and modify them as necessary.

To reproduce the manuscript results, use the same input data, reference sequence, software environment, coordinate definitions, and analysis settings described in the manuscript.

## Main output files

Depending on the platform and selected analysis, the scripts generate files including:

```text
fig/QT/sig_compare_all.png
fig/QT/sig_compare_font.png
fig/QT/sig_compare_dAAL.png
fig/QT/sig_distribution_font.png
fig/QT/sig_distribution_dAAL.png
fig/QT/zscore_threshold_1_out_ratio.png
fig/QT/zscore_threshold_2_out_ratio.png
fig/QT/zscore_threshold_3_out_ratio.png
fig/QT/zscore_threshold_1_confusion_matrix.png
fig/QT/zscore_threshold_2_confusion_matrix.png
fig/QT/zscore_threshold_3_confusion_matrix.png
fig/QT/roc_curve.png
fig/QT/hotmap.png
fig/QT/DIS_radar_plot.png
```

Equivalent signal-comparison and Z-score results are generated under:

```text
fig/ONT/
```

The sequencing-error analysis generates files including:

```text
data/ONT_odds_ratio.csv
fig/ONT/ONT_odds_ratio.png
fig/ONT/ONT_esb_radar.png
fig/ONT/ONT_quality.png
```

When the platform variable is changed to QT, corresponding output files are generated for the QT dataset.

The exact output files depend on the platform and configuration selected inside each script.

## Data availability

The raw ionic-current signal files, basecalled and reference-aligned sequencing reads, BAM index files, prepared reference-aligned signal matrices, custom reference sequence, and processed analysis results are publicly available in Zenodo:

https://doi.org/10.5281/zenodo.21740475

The dataset includes data generated using both Oxford Nanopore and Qitan Technology sequencing platforms.

Large research-data files are archived in Zenodo and are therefore not stored directly in this GitHub repository.

## Code availability

All custom downstream scripts used for nanopore signal comparison, Z-score-based detection, statistical evaluation, sequencing-error analysis, classification-performance evaluation, and figure generation are publicly available in this GitHub repository:

https://github.com/liuhaha0607/nanopore-zscore-analysis

Version 1.0.0 of the repository has been permanently archived in Zenodo:

https://doi.org/10.5281/zenodo.21737768

The DOI badge at the top of this README points to the latest archived release of the repository.

## Citation

### Software

If you use the analysis code, please cite the archived software version:

> Liu, Ran. *Nanopore Signal Analysis and Z-score-based Detection* (Version 1.0.0) [Software]. Zenodo. https://doi.org/10.5281/zenodo.21737768

### Dataset

If you use the associated data, please cite:

> Liu, Ran, Shen, Yezhuang, Mao, Jie, Li, Chunzheng, Zhang, Yawei, Hu, Mandong, Zhang, Yizhe, Han, Qiuying, Gong, Weili, Chen, Liang, He, Kun, Zhou, Tao, Li, Weihua, and Xie, Xianxing. *Nanopore signal and sequencing data for control and dA-AL-II-modified synthetic DNA generated using Oxford Nanopore and Qitan Technology platforms* (Version 1.0.0) [Dataset]. Zenodo. https://doi.org/10.5281/zenodo.21740475

The version-specific DOI identifies the exact archived software or dataset release used for the study.

## License

This project is distributed under the MIT License.

See the `LICENSE` file for the full license text.

The MIT License applies to the original code contained in this repository.
