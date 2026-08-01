# Nanopore Signal Analysis and Z-score-based Detection

## Overview

This repository contains the custom downstream analysis scripts used for nanopore signal comparison, Z-score-based detection, sequencing-error analysis, statistical testing, classification-performance evaluation, and figure generation.

The repository includes scripts for:

* comparison and visualization of reference-aligned nanopore signal levels;
* position-specific statistical analysis of signal differences;
* Z-score-based classification of control and dA-AL-II reads;
* evaluation of different Z-score thresholds;
* confusion-matrix and receiver operating characteristic (ROC) analysis;
* simulated mixture-proportion analysis;
* sequencing-error, odds-ratio, and Fisher's exact-test analysis;
* generation of the figures and tabular results used in the study.

The signal-analysis scripts begin with prepared, reference-aligned signal matrices. This repository does not provide an end-to-end implementation for regenerating all aligned signal matrices directly from raw nanopore signal files.

---

## Repository contents

### `plot_signal_compare.py`

Compares reference-aligned signal traces between control and dA-AL-II samples.

The script generates:

* whole-region signal-comparison plots;
* local signal-comparison plots;
* position-specific kernel-density plots;
* separate results for QT and ONT datasets.

A median filter is used for visualization of the aligned signal traces.

### `zscore_read_classifier.py`

Performs the position-specific Z-score-based detection analysis.

The script:

* uses part of the control reads to calculate the background mean and standard deviation at each reference position;
* calculates the absolute Z-score of each control and dA-AL-II read;
* evaluates Z-score thresholds of 1, 2, and 3;
* calculates false-positive rate, true-positive rate, precision, and F1 score;
* generates normalized confusion matrices;
* generates an ROC curve and calculates the area under the curve;
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

### `generate_odds.py`

Performs position-specific sequencing-error analysis using aligned BAM files and a reference FASTA file.

The script:

* calculates mismatch, insertion, and deletion counts;
* calculates total sequencing-error rates;
* summarizes base-quality scores;
* compares the dA-AL-II and control samples;
* calculates odds ratios using a pseudocount;
* performs one-sided Fisher's exact tests;
* generates odds-ratio, error-rate radar, and quality-comparison plots;
* exports the position-specific statistical results as a CSV file.

---

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

The figure directories are created automatically by the scripts when the analyses are run.

---

## Reproducibility scope and upstream signal processing

Raw nanopore signals were processed, normalized, and aligned to reference positions before the custom downstream analyses provided in this repository.

Upstream signal processing used established nanopore-processing software and platform-specific resources. For the QT platform, reference-level signal mapping depended on a proprietary theoretical-current maptable supplied by Qitan Technology.

The proprietary QT maptable and its associated upstream implementation are subject to third-party restrictions. The authors do not have permission to publicly redistribute these materials. They are therefore not included in this repository.

Consequently, this repository supports reproduction of the custom downstream analyses from compatible, prepared reference-aligned signal matrices. It does not support complete regeneration of the QT aligned signal matrices directly from the raw QT signal files.

The proprietary restriction applies only to the upstream QT signal-mapping resource. All custom downstream scripts developed for signal comparison, statistical testing, Z-score-based detection, classification evaluation, sequencing-error analysis, and figure generation are publicly provided in this repository.

---

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

Install the required Python packages using:

```bash
pip install -r requirements.txt
```

Alternatively, create the provided Conda environment:

```bash
conda env create -f environment.yml
conda activate nanopore-zscore-analysis
```

Using the provided environment specification is recommended for reproducing the computational results.

---

## Preparing the input directory

Create a directory named `data` in the repository root before running the scripts:

```bash
mkdir data
```

The expected local directory structure is:

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

Only files required for the analyses being performed need to be placed in the directory.

---

## Expected signal-matrix inputs

The following prepared signal-matrix files are used by the signal-comparison, DIS, and Z-score scripts:

```text
data/QT_control_align_res.pkl
data/QT_dAAL_align_res.pkl
data/ONT_control_align_res.pkl
data/ONT_dAAL_align_res.pkl
```

Each pickle file is expected to contain a Python dictionary with a key named:

```python
"sigs"
```

The value associated with `sigs` must be a two-dimensional signal matrix in which:

* each row represents an individual nanopore read;
* each column represents a reference-aligned signal measurement;
* the current downstream scripts use six signal measurements per reference position.

The compatible signal matrices must have already undergone the required upstream signal processing, normalization, reference alignment, and resampling.

The proprietary QT theoretical-current maptable and the upstream procedure used to generate the QT aligned signal matrices are not distributed in this repository.

---

## Expected sequencing-analysis inputs

The sequencing-error analysis expects:

```text
data/ONT_control.bam
data/ONT_dAAL.bam
data/ref.fasta
```

Indexed BAM files are recommended:

```text
data/ONT_control.bam.bai
data/ONT_dAAL.bam.bai
```

The BAM files must be aligned to the reference sequence contained in `ref.fasta`.

The reference sequence and reference names used in the BAM and FASTA files must be consistent.

---

## Running the analyses

Run the scripts from the root directory of the repository.

### Signal comparison and visualization

```bash
python plot_signal_compare.py
```

This script processes the QT and ONT aligned signal matrices and saves figures under:

```text
fig/QT/
fig/ONT/
```

### Z-score-based detection

```bash
python zscore_read_classifier.py
```

This script processes the QT and ONT aligned signal matrices and generates:

* Z-score threshold plots;
* confusion matrices;
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

The platform used in this analysis is selected inside the script.

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

The script reads the aligned BAM files and reference FASTA file and generates:

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

---

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
* target position;
* statistical significance threshold;
* output directories.

The current signal-analysis scripts use six aligned signal measurements per reference position.

The main highlighted target position is reference position 506.

Before running an analysis on another dataset, review the configuration variables in the corresponding script and update them as necessary.

To reproduce the manuscript results, use the same input data, reference sequence, software environment, coordinate definitions, and analysis settings described in the manuscript.

---

## Output files

Depending on the selected platform and analysis, the scripts generate files including:

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

Equivalent signal and Z-score results may be generated under:

```text
fig/ONT/
```

The sequencing-error analysis additionally generates files including:

```text
data/ONT_odds_ratio.csv
fig/ONT/ONT_odds_ratio.png
fig/ONT/ONT_esb_radar.png
fig/ONT/ONT_quality.png
```

The exact outputs depend on the platform and configuration selected inside each script.

---

## Important limitations

This repository does not include:

* the proprietary QT theoretical-current maptable;
* proprietary QT signal-mapping software or implementation details;
* raw `.fast5`, `.pod5`, or `.h5` signal files;
* large BAM files;
* large processed signal matrices;
* confidential or third-party materials that the authors are not authorized to redistribute.

The absence of the proprietary QT maptable means that the QT aligned signal matrices cannot be regenerated from raw QT signal files using this repository alone.

The public scripts nevertheless document and implement the custom downstream computational procedures used for:

* signal visualization;
* signal-difference testing;
* position-specific Z-score calculation;
* threshold-based classification;
* ROC and confusion-matrix analysis;
* mixture-proportion evaluation;
* sequencing-error analysis;
* odds-ratio calculation;
* figure generation.

---

## Data availability

The raw nanopore signal files and associated sequencing data generated in this study will be deposited in a public sequencing archive before publication.

Replace the placeholders below after the archive submission has been approved:

```text
Archive: [NCBI Sequence Read Archive or European Nucleotide Archive]
BioProject/Study accession: [ADD BIOPROJECT OR STUDY ACCESSION]
Run accession(s): [ADD RUN ACCESSION NUMBERS]
```

The proprietary QT theoretical-current maptable is owned or controlled by a third party and cannot be redistributed by the authors.

Raw signal files, BAM files, and large processed signal matrices are not stored directly in this GitHub repository.

---

## Code availability

All custom downstream scripts used for nanopore signal comparison, Z-score-based detection, statistical evaluation, sequencing-error analysis, classification-performance evaluation, and figure generation are publicly available at:

```text
https://github.com/liuhaha0607/nanopore-zscore-analysis
```

A versioned release of this repository will be archived in Zenodo.

After the Zenodo record has been created, replace the placeholder below:

```text
Zenodo DOI: [ADD ZENODO DOI]
```

The public availability of the custom downstream code does not extend to the proprietary QT theoretical-current maptable or associated third-party upstream implementation.

---

## Citation

Until a Zenodo DOI is available, this repository may be cited as:

```text
Liu, Ran. Nanopore Signal Analysis and Z-score-based Detection.
GitHub repository:
https://github.com/liuhaha0607/nanopore-zscore-analysis
```

After creating the Zenodo release, update this section with the version number and DOI:

```text
Liu, Ran. Nanopore Signal Analysis and Z-score-based Detection,
version 1.0.0. Zenodo. https://doi.org/[ADD ZENODO DOI]
```

---

## License

This project is distributed under the MIT License.

See the `LICENSE` file for the full license text.

The MIT License applies only to the original code contained in this repository. It does not grant rights to any proprietary third-party model, maptable, software, data, or other restricted material referred to in the documentation.
