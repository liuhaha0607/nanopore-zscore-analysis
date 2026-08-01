# Nanopore Signal Analysis and Z-score-based Detection

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

The signal-analysis scripts begin with prepared, normalized, reference-aligned signal matrices. This repository does not provide an end-to-end implementation for regenerating these matrices directly from raw nanopore signal files.

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
* calculates combined mismatch and deletion error rates;
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

The `data/` directory is not included in the public repository because it may contain large raw, aligned, or processed research files.

The `fig/` subdirectories are created automatically when the scripts are run.

## Reproducibility scope and upstream signal processing

Raw nanopore signals were processed, normalized, and aligned to reference positions before the custom downstream analyses provided in this repository.

Upstream signal processing was performed outside this repository using platform-specific workflows and resources. For the QT platform, reference-level signal mapping depended on a proprietary theoretical-current maptable supplied by Qitan Technology.

The proprietary QT maptable and its associated upstream implementation are subject to third-party restrictions. The authors do not have permission to publicly redistribute these materials, and they are therefore not included in this repository.

Consequently, this repository supports reproduction of the custom downstream analyses from compatible, prepared reference-aligned signal matrices. It does not support complete regeneration of the QT aligned signal matrices directly from raw QT signal files.

The proprietary restriction applies only to the upstream QT signal-mapping resource. All custom downstream scripts developed for signal comparison, statistical testing, Z-score-based detection, classification evaluation, sequencing-error analysis, and figure generation are publicly provided in this repository.

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

## Preparing the input directory

Create a directory named `data` in the repository root before running the scripts:

```bash
mkdir data
```

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

The proprietary QT theoretical-current maptable and the upstream procedure used to generate the QT aligned signal matrices are not distributed in this repository.

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

The current script uses the ONT BAM files and generates:

* a position-specific odds-ratio table;
* an odds-ratio plot;
* an error-rate radar plot;
* a quality-score comparison plot.

The CSV result is saved under:

```text
data/
```

The figures are saved under:

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

Equivalent signal-comparison and Z-score results are generated under `fig/ONT/`.

The sequencing-error analysis generates:

```text
data/ONT_odds_ratio.csv
fig/ONT/ONT_odds_ratio.png
fig/ONT/ONT_esb_radar.png
fig/ONT/ONT_quality.png
```

The exact output files depend on the platform and configuration selected inside each script.

## Important limitations

This repository does not include:

* the proprietary QT theoretical-current maptable;
* proprietary QT signal-mapping software or implementation details;
* raw `.fast5`, `.pod5`, or `.h5` signal files;
* large BAM files;
* large processed signal matrices;
* confidential or third-party materials that the authors are not authorized to redistribute.

The absence of the proprietary QT maptable means that the QT aligned signal matrices cannot be regenerated from raw QT signal files using this repository alone.

The publicly available scripts nevertheless document and implement the custom downstream computational procedures used for:

* signal visualization;
* signal-difference testing;
* position-specific Z-score calculation;
* threshold-based classification;
* ROC and confusion-matrix analysis;
* mixture-proportion evaluation;
* mismatch and deletion error analysis;
* odds-ratio calculation;
* figure generation.

## Data availability

The raw nanopore signal files and associated sequencing data generated in this study will be deposited in a public sequencing archive before publication.

Accession information will be added after the archive submission has been processed:

```text
Archive: [NCBI Sequence Read Archive or European Nucleotide Archive]
BioProject/Study accession: [ADD BIOPROJECT OR STUDY ACCESSION]
Run accession(s): [ADD RUN ACCESSION NUMBERS]
```

The proprietary QT theoretical-current maptable is owned or controlled by a third party and cannot be redistributed by the authors.

Raw signal files, BAM files, and large processed signal matrices are not stored directly in this GitHub repository.

## Code availability

All custom downstream scripts used for nanopore signal comparison, Z-score-based detection, statistical evaluation, sequencing-error analysis, classification-performance evaluation, and figure generation are publicly available at:

https://github.com/liuhaha0607/nanopore-zscore-analysis

A versioned release of this repository will be archived in Zenodo.

After the Zenodo record has been created, the DOI will be added here:

```text
Zenodo DOI: [ADD ZENODO DOI]
```

The public availability of the custom downstream code does not extend to the proprietary QT theoretical-current maptable or associated third-party upstream implementation.

## Citation

Until a Zenodo DOI is available, cite this repository as:

```text
Liu, Ran. Nanopore Signal Analysis and Z-score-based Detection.
GitHub repository: https://github.com/liuhaha0607/nanopore-zscore-analysis
```

After the Zenodo release has been created, update the citation with the version number and DOI:

```text
Liu, Ran. Nanopore Signal Analysis and Z-score-based Detection,
version 1.0.0. Zenodo. https://doi.org/[ADD ZENODO DOI]
```

## License

This project is distributed under the MIT License.

See the `LICENSE` file for the full license text.

The MIT License applies only to the original code contained in this repository. It does not grant rights to any proprietary third-party model, maptable, software, data, or other restricted material referred to in this documentation.
