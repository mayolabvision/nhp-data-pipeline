# nhp-data-pipeline

Custom neural and behavioral data processing pipeline for the **Mayo/Smith Lab** at the University of Pittsburgh. Designed to handle multi-probe recordings from non-human primates, automate spike sorting and postprocessing, and compile analysis-ready tables — all running on the **H2P CRC Pitt cluster** via SLURM.

---

## Overview

This pipeline takes raw electrophysiology and behavior recordings from the lab's NHP setups and processes them through a standardized, reproducible workflow. Processing parameters are defined in JSON **protocol files**, and outputs are organized using content-based hashes so that results are always traceable back to the exact parameters that produced them.

The pipeline handles two classes of neural recordings:

- **Neuropixel probes** — recorded with IMEC hardware via SpikeGLX (`.bin` / `.meta` format)
- **Plexon / FHC probes** — recorded with a Ripple Grapevine system (`.ns5` / `.nev` format)

Behavioral data (eye position, task events, reward signals, etc.) is also recorded with the **Ripple Grapevine** and processed alongside the neural data.

---

## Pipeline Stages

Each stage is submitted as a separate SLURM job, allowing flexible resource allocation (CPUs for preprocessing, GPUs for sorting, etc.).

### 1. Preprocessing — `submit_preprocess.sh`
Loads raw data for a given session and probe, converts Ripple files to binary if needed, and applies a preprocessing chain via [SpikeInterface](https://spikeinterface.readthedocs.io/en/latest/). Steps include bandpass filtering, common-median referencing, and motion correction. Outputs are saved under a hash derived from the preprocessing parameters.

> Runs on: **SMP cluster, high-mem partition** (16 CPUs)

### 2. Shake Trimming — `submit_trimming.sh`
Estimates probe drift across the recording and flags or removes epochs where motion exceeds a stability threshold. This is an optional step controlled by the protocol file.

> Runs on: **SMP cluster**

### 3. Spike Sorting — `submit_sorting.sh`
Runs **Kilosort 4** (via SpikeInterface) on the preprocessed recording. GPU-accelerated. Supports per-probe parallelism via SLURM array jobs.

> Runs on: **GPU cluster, A100 partition** (2 GPUs, 12 CPUs)

### 4. Postprocessing — `submit_postprocess.sh`
Builds a SpikeInterface `SortingAnalyzer` from the sorting output. Computes waveforms, templates, principal components, correlograms, and unit quality metrics. Saves a protocol-to-hash mapping for downstream tracking.

> Runs on: **SMP cluster** (8 CPUs)

### 5. Visualization — `submit_widgets.sh` / `submit_plots.sh`
Generates summary figures and SpikeInterface widgets for quality review, including probe motion traces, noise levels, preprocessing steps, and per-unit plots.

### 6. MATLAB Table Compilation — `submit_matlabTbl.sh`
Final stage. Runs `PROCESS_RECORDING.m` to parse sorting output and behavior data into a unified `.mat` table per session. Handles spike time binning, event alignment, LFP extraction, sync pulse matching (Ripple ↔ Neuropixel), and task-specific formatting.

> Runs on: **SMP cluster, high-mem partition** (32 CPUs, MATLAB R2023a)

---

## Recording Profiles

The pipeline uses an abstract `RecordingProfile` class with hardware-specific subclasses:

| Profile | Hardware | File Format |
|---|---|---|
| `NeuropixelProfile` | IMEC Neuropixel | SpikeGLX `.bin` / `.meta` |
| `PlexonProfile` | Ripple Grapevine + Plexon headstage | `.ns5` / `.nev` |
| `FHCProfile` | Ripple Grapevine + FHC microdrive | `.ns5` / `.nev` |
| `BehaviorProfile` | Ripple Grapevine (behavior only) | `.ns5` / `.nev` |

The correct profile is selected automatically based on session metadata.

---

## Protocol Files

Processing parameters are defined in `.json` protocol files stored at `/ix1/pmayo/protocols/`. A protocol specifies the preprocessing chain, motion correction settings, spike sorter parameters, and optional trimming thresholds. Example protocols live in `examples/protocols/`.

Outputs are organized by a hash computed from the protocol parameters, ensuring that results are reproducible and parameter changes automatically produce new output directories.

---

## Repository Structure

```
nhp-data-pipeline/
├── main_pipeline.py          # Entry points: preprocess, sort, postprocess, etc.
├── config.py                 # Cluster paths (data, envs, packages)
├── pipeline/
│   ├── base.py               # Abstract RecordingProfile base class
│   ├── neuropixel.py         # Neuropixel-specific processing
│   ├── plexon.py             # Plexon/Ripple processing
│   ├── fhc.py                # FHC/Ripple processing
│   ├── behavior.py           # Behavior-only processing
│   ├── si_tools.py           # SpikeInterface utilities
│   ├── si_plots.py           # SpikeInterface visualization
│   ├── ks_tools.py           # Kilosort output utilities
│   ├── path_utils.py         # Hash-based path management
│   └── catgt_utils.py        # CatGT preprocessing utilities
├── matlab/
│   ├── PROCESS_RECORDING.m   # Main MATLAB entry point
│   ├── parse_SortingToTbl.m  # Parse sorting output into table
│   ├── extract_lfpData.m     # LFP extraction
│   ├── match_syncPulses_RipToNP.m  # Ripple ↔ Neuropixel sync
│   └── ...                   # Additional helper scripts
├── submit_preprocess.sh      # SLURM: preprocessing
├── submit_trimming.sh        # SLURM: shake trimming
├── submit_sorting.sh         # SLURM: spike sorting (GPU)
├── submit_postprocess.sh     # SLURM: postprocessing + metrics
├── submit_widgets.sh         # SLURM: SpikeInterface widgets
├── submit_plots.sh           # SLURM: summary figures
├── submit_matlabTbl.sh       # SLURM: MATLAB table compilation
├── examples/
│   ├── protocols/            # Example protocol JSON files
│   ├── metadata/             # Example session metadata files
│   └── probes/               # Example probe geometry files
└── notebooks/                # Development and testing notebooks
```

---

## Cluster Setup

All jobs run on the **University of Pittsburgh H2P CRC cluster**. Data lives at `/ix1/pmayo/lab_NHPdata/`. The conda environment is at `/ix1/pmayo/envs/NHPipe`.

To submit a full pipeline run for a session:

```bash
# 1. Preprocess (per probe, array job)
sbatch submit_preprocess.sh <SESSION> <PROTOCOL>

# 2. Spike sort (GPU, per probe)
sbatch submit_sorting.sh <SESSION> <PROTOCOL>

# 3. Postprocess + quality metrics
sbatch submit_postprocess.sh <SESSION> <PROTOCOL>

# 4. Compile MATLAB table
sbatch submit_matlabTbl.sh <SESSION> <PROTOCOL>
```

---

## Coming Soon

- Automated pipeline chaining (submit full workflow as a DAG with SLURM dependencies)
- Drift-corrected LFP extraction for Neuropixel recordings
- More thorough debugging jupyter notebook
- Full installation/debugging tutorials
- "Quick" versions of the pipeline for behavior only 
- Improved probe geometry tooling for custom FHC array configurations

---

## Questions?

For questions about this pipeline, contact **Kendra Noneman** — knoneman@andrew.cmu.edu
