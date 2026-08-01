# Nanopore Signal Analysis and Z-score-based Detection

This repository contains the analysis scripts used for nanopore signal comparison,
Z-score-based detection, sequencing-error analysis, statistical testing, and figure
generation.

## Scripts

- `dis_plot.py`  
  Calculates site-level signal differences between treatment and control samples,
  performs Welch's t-tests, calculates standardized signal deviations, and generates
  a radar plot.

- `generate_odds.py`  
  Calculates position-specific sequencing error rates from BAM files, compares test
  and control samples using odds ratios and Fisher's exact tests, summarizes quality
  scores, and generates odds-ratio, radar, and quality-comparison plots.

- `plot_signal_compare.py`  
  Compares aligned signal traces from control and dA-AL-II samples and generates
  whole-region, local-region, and kernel-density signal plots for QT and ONT data.

- `zscore_read_classifier.py`  
  Calculates absolute Z-scores using the control signal distribution, evaluates
  multiple detection thresholds, generates confusion matrices and ROC curves, and
  evaluates predicted positive proportions in simulated sample mixtures.

## Repository structure

```text
nanopore-zscore-analysis/
├── README.md
├── requirements.txt
├── environment.yml
├── .gitignore
├── dis_plot.py
├── generate_odds.py
├── plot_signal_compare.py
├── zscore_read_classifier.py
├── data/
│   └── README.md
└── fig/
```

The raw and processed research data are not stored in this GitHub repository.
Raw nanopore data should be deposited in a public sequencing archive and cited
using its accession number.

## Software requirements

A recent Python 3 installation is recommended. Install the required packages with:

```bash
pip install -r requirements.txt
```

Alternatively, create a Conda environment with:

```bash
conda env create -f environment.yml
conda activate nanopore-zscore-analysis
```

## Expected input files

Place the required local analysis files in the `data/` directory before running
the scripts. See `data/README.md` for the expected filenames.

## Running the analysis

From the repository directory, run the required scripts individually:

```bash
python plot_signal_compare.py
python zscore_read_classifier.py
python dis_plot.py
python generate_odds.py
```

The scripts save figures under `fig/QT/` or `fig/ONT/` and save tabular results
under `data/`, depending on the script.

## Important configuration notes

Several analysis settings are currently defined inside the scripts, including:

- platform (`QT` or `ONT`);
- signal resampling factor;
- analyzed sequence positions;
- Z-score thresholds;
- minimum read depth;
- highlighted target position;
- input filenames.

To reproduce the manuscript results, use the same settings reported in the
Methods section and the same software environment described in this repository.

## Data availability

Raw nanopore signal files and associated sequencing data are available from:

- Archive: NCBI Sequence Read Archive (SRA) or European Nucleotide Archive (ENA)
- BioProject/Study accession: `[ADD ACCESSION NUMBER]`
- Run accession(s): `[ADD RUN ACCESSION NUMBERS]`

Do not upload raw `.fast5`, `.pod5`, `.h5`, `.bam`, or large processed-data files
directly to this GitHub repository.

## Code citation

A versioned release of this repository should be archived in Zenodo. After Zenodo
creates a DOI, replace the placeholder below:

> Code DOI: `[ADD ZENODO DOI]`

## License

Before public release, add an open-source license approved by all authors and
your institution. The MIT License is commonly used for research software, but
the appropriate license should be confirmed by the authors.
