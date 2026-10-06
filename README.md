# Wideband RF/Photonic Beamforming Simulator

System-level beamforming project implemented in MATLAB and Ansys Lumerical INTERCONNECT.

The project compares fixed phase-shifter beamforming, ideal true-time-delay (TTD), and photonic TTD for an 8-element ULA centered at 10 GHz.

## Contents

- `MATLAB/`  
  Array-factor simulations, beam squint, spatial sampling, amplitude weighting, multibeam operation, RF non-idealities, wideband PS vs TTD sweeps, and end-to-end comparison.

- `INTERCONNECT/`  
  8-channel photonic TTD model, ideal/non-ideal/thermal setup scripts, and export script.

- `Beamforming Project Report.pdf`  
  Full project report with methodology, figures, limitations, and final quantitative comparison.

## Main Results

For a 30° steering target over approximately 7.5–12.5 GHz:

- Phase shifter maximum pointing error: 11.47°
- Non-ideal photonic TTD maximum pointing error: 0.59°
- Thermal + noise photonic TTD maximum pointing error: 1.31°

## Requirements

- MATLAB
- Ansys Lumerical INTERCONNECT

The exported `.mat` files are included so the MATLAB end-to-end analysis can be run without rerunning INTERCONNECT.

## Usage

Run the MATLAB scripts individually for the RF analyses.

For the end-to-end comparison, run:

`MATLAB/end_to_end_comparison.m`

For INTERCONNECT, open:

`INTERCONNECT/beamforming_project.icp`

and use the corresponding `.lsf` setup and export scripts.
