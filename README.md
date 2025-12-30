# SAMS - Semi-Automatic MEA Spike Sorting

A MATLAB-based pipeline for semi-automatic spike sorting and network burst analysis of multi-electrode array (MEA) recordings.

> **Note:** SAMS is currently developed and tested for **Windows only**. macOS/Linux support may be added in future versions.

## Features

- **Automated Spike Sorting**: Template matching and clustering algorithms for spike classification
- **Network Burst Detection**: Identification and analysis of network-level bursting activity
- **Individual Unit Analysis**: Comprehensive metrics for sorted units
- **Interactive GUI**: User-friendly interface for manual curation and visualization
- **Batch Processing**: Process multiple wells and electrodes efficiently
- **Export Capabilities**: Generate reports and export results to PowerPoint

## Requirements

- **Windows OS** (tested on Windows 10/11)
- [MATLAB Runtime R2023a](https://www.mathworks.com/products/compiler/mcr/index.html)

## Installation

**Automatic Sorting:**
1. Download and install [MATLAB Runtime R2023a](https://www.mathworks.com/products/compiler/mcr/index.html) if not already installed
2. Download the **entire** `SemiAutomaticMEASpikeSorting/for_redistribution_files_only/` folder (all files are required, including `SAMS1.pptx` for output slides and `splash.png`)
3. Run `SemiAutomaticMEASpikeSorting.exe`

**Manual Curation (requires output from Automatic Sorting):**
1. Download the **entire** `manualspikesorting/for_redistribution_files_only/` folder (all files are required)
2. Run `manualspikesorting.exe`
3. Load the output files from automatic sorting to review and curate results

## Usage

### Quick Start

1. Launch the automatic sorting application (GUI)
2. Select folder that contains `.spk` files (Axion BioSystems format)
3. (Optional) Select wells to process (default: all wells)
4. Run automatic sorting
5. Review (and manually curate results in manual sorting app as needed)
6. Export results

See [SAMS_User_Manual.docx](docs/SAMS_User_Manual.docx) for details.

## Project Structure

```
SAMS/
├── README.md
├── LICENSE
├── THIRD_PARTY_NOTICES.md
├── SemiAutomaticMEASpikeSorting/    # Automatic sorting standalone app
│   └── for_redistribution_files_only/
├── manualspikesorting/              # Manual curation standalone app
│   └── for_redistribution_files_only/
├── src/                             # Shared source code
│   ├── automatic_sorting.mlapp     # Automatic sorting GUI
│   ├── manual_sorting_03252025.mlapp  # Manual sorting GUI
│   ├── AxionFileLoader/            # Axion file reading library
│   └── [additional source files]
├── sample files/                    # Sample data for testing
│   ├── Axion Spike File/           # Sample .spk file
│   ├── SAMS Output Files/          # SAMS results for the sample
│   ├── For OFS - Unsorted/         # Unsorted .nex/.plx files
│   └── For OFS - Sorted/           # Sorted .nex/.plx files
└── docs/                            # Documentation
    └── SAMS_User_Manual.docx
```

## Sample Files

The `sample files/` folder contains test data with 3 wells of FXS2 neurons (24 active electrodes) and 3 wells of CTRL2 neurons (16 active electrodes):

- **Axion Spike File/** - Sample `.spk` file from Axion BioSystems
- **SAMS Output Files/** - SAMS processing results for the sample
- **For OFS - Unsorted/** - Unsorted `.nex` and `.plx` files (converted using Axion's Data Export Tool)
- **For OFS - Sorted/** - Manually sorted files using Offline Sorter (OFS)

## Documentation

- [SAMS_User_Manual.docx](docs/SAMS_User_Manual.docx) - Comprehensive user guide

## Third-Party Components

This project uses the following third-party components:
- **AxionFileLoader** - MIT License, (c) Axion BioSystems
- **InterX** - Curve intersection function by NS
- **Hartigan's Dip Test** - Algorithm AS 217 implementation by F. Mechler

See [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) for full details.

## Authors

Xiaoxuan Ren, Carissa L. Sirois, Raymond Doudlah, Ethan Dayley, Natasha M. Mendez-Albelo, Aviad Hai, Ari Rosenberg, Xinyu Zhao

## License

See the [LICENSE](LICENSE) file for details.

## Citation

If you use this software in your research, please cite:

> Ren X, Sirois CL, Doudlah R, Mendez-Albelo NM, Hai A, Rosenberg A, Zhao X. (2025). A Semi-Automated MEA Spike sorting (SAMS) method for high throughput assessment of cultured neurons. *bioRxiv*. DOI: [10.1101/2025.02.08.637245](https://doi.org/10.1101/2025.02.08.637245)

## Support

For questions or issues, please contact xinyu.zhao@wisc.edu
