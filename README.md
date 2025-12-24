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

### Option 1: Standalone Application (Recommended for end users)

1. Download the latest release from the `SemiAutomaticMEASpikeSorting/for_redistribution_files_only/` folder
2. Run `SemiAutomaticMEASpikeSorting.exe`
3. If prompted, install the MATLAB Runtime (will be downloaded automatically)

### Option 2: Running from Source

1. Clone this repository:
   ```
   git clone https://github.com/YOUR_USERNAME/SAMS.git
   ```
2. Add the `src/` folder and subfolders to your MATLAB path:
   ```matlab
   addpath(genpath('path/to/SAMS/src'));
   ```
3. Run the main application:
   ```matlab
   run_SAMS
   ```
   Or launch the GUI:
   ```matlab
   automatic_sorting
   ```

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
├── SemiAutomaticMEASpikeSorting/    # Compiled standalone application
│   ├── for_redistribution/          # Web installer
│   ├── for_redistribution_files_only/  # Standalone exe
│   └── for_testing/                 # Development build
├── src/                             # Source code
│   ├── run_SAMS.m                   # Main entry point
│   ├── process_file.m              # File processing
│   ├── process_electrode.m         # Electrode processing
│   ├── automatic_sorting.mlapp     # Main GUI
│   ├── manual_sorting_03252025.mlapp  # Manual sorting GUI
│   ├── AxionFileLoader/            # Axion file reading library
│   └── [additional source files]
└── docs/                            # Documentation
    ├── SAMS_User_Manual.docx
    └── SAMS1.pptx
```

## Documentation

See the `docs/` folder for:
- `SAMS_User_Manual.docx` - Comprehensive user guide
- `SAMS1.pptx` - Overview presentation

## Third-Party Components

This project uses the following third-party components:
- **AxionFileLoader** - MIT License, (c) Axion BioSystems
- **InterX** - Curve intersection function by NS
- **Hartigan's Dip Test** - Algorithm AS 217 implementation by F. Mechler

See [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) for full details.

## Author

Xiaoxuan Ren - Zhao Lab

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this software in your research, please cite:

```
[Citation information to be added]
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Support

For questions or issues, please open an issue on GitHub.
