# SAMS - Semi-Automatic MEA Spike Sorting

A MATLAB-based pipeline for semi-automatic spike sorting and network burst analysis of multi-electrode array (MEA) recordings.

## Features

- **Automated Spike Sorting**: Template matching and clustering algorithms for spike classification
- **Network Burst Detection**: Identification and analysis of network-level bursting activity
- **Individual Unit Analysis**: Comprehensive metrics for sorted units
- **Interactive GUI**: User-friendly interface for manual curation and visualization
- **Batch Processing**: Process multiple wells and electrodes efficiently
- **Export Capabilities**: Generate reports and export results to PowerPoint

## Requirements

- MATLAB R2023a or later
- Statistics and Machine Learning Toolbox
- Signal Processing Toolbox

## Installation

**Automatic Sorting:**
1. Download from `SemiAutomaticMEASpikeSorting/for_redistribution_files_only/`
2. Run `SemiAutomaticMEASpikeSorting.exe`
3. If prompted, install the MATLAB Runtime (will be downloaded automatically)

**Manual Curation (requires output from Automatic Sorting):**
1. Download from `manualspikesorting/for_redistribution_files_only/`
2. Run `manualspikesorting.exe`
3. Load the output files from automatic sorting to review and curate results

## Usage

### Quick Start

1. Launch the application (GUI or command line)
2. Load your `.spk` file (Axion BioSystems format)
3. Select wells and electrodes to process
4. Run automatic sorting
5. Review and manually curate results as needed
6. Export results

### Input Data

This pipeline is designed to work with spike files (`.spk`) from Axion BioSystems MEA platforms.

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
│   ├── run_SAMS.m                   # Main entry point
│   ├── automatic_sorting.mlapp     # Automatic sorting GUI
│   ├── manual_sorting_03252025.mlapp  # Manual sorting GUI
│   ├── AxionFileLoader/            # Axion file reading library
│   └── [additional source files]
└── docs/                            # Documentation
    └── SAMS_User_Manual.docx
```

## Documentation

See the `docs/` folder for:
- `SAMS_User_Manual.docx` - Comprehensive user guide

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
