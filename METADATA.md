# METADATA — Dataset Descriptions and Variable Definitions

## Overview

This dataset comprises two MATLAB scripts (`demo1.m` and `demo2.m`) that implement and evaluate signal denoising algorithms, with a focus on the **Recursive Sliding-window Band-Adaptive Variational Mode Decomposition (RSBAVMD)** method. The scripts process both **synthetic (simulated)** signals and **real measured** vibration data, and compare RSBAVMD against classical VMD, EMD, and Butterworth bandpass filtering.

- **Associated data file**: `data.mat` — contains a field `Sig`, which is a real measured vibration signal time series.
- **Script language**: MATLAB (tested with Signal Processing Toolbox; requires `emd` from a compatible toolbox).
- **Research domain**: Structural health monitoring / vibration signal denoising / time-frequency analysis.

---

## 1. File: `demo1.m` — Simulated & Measured Signal Denoising Comparison

### 1.1 Experiment Purpose

`demo1.m` performs a comprehensive denoising benchmark on a **synthetic multi-component signal**.

### 1.2 Simulated Signal Variables

| Variable     | Type       | Unit  | Definition                                                                |
|--------------|------------|-------|---------------------------------------------------------------------------|
| `fs`         | scalar     | Hz    | Sampling frequency. Fixed at **50 Hz**.                                   |
| `t`          | [N x 1]    | s     | Equally-spaced time vector from 0 to T-1/fs in steps of 1/fs.             |
| `s_d`        | [N x 1]    | —     | **Damped vibration mode**                                                 |
| `s_1`        | [N x 1]    | —     | **Intermodulation component 1**                                           |
| `s_2`        | [N x 1]    | —     | **Intermodulation component 2**                                           |
| `x_clean`    | [N x 1]    | —     | **Composite clean signal**  Ground truth for denoising evaluation.        |
| `snr_input`  | scalar     | dB    | **Input SNR level** for noise addition. Set to **2 dB**.                  |
| `x_noisy`    | [N x 1]    | —     | **Noisy signal** . White Gaussian noise.                                  |

### 1.3 RSBAVMD Parameters (`demo1.m`, simulated part)

| Variable         | Type    | Unit | Definition                                               |
|------------------|---------|------|----------------------------------------------------------|
| `Sig`            | [N x 1] | —    | Input signal to RSBAVMD.                                 |
| `SigFFT`         | [N x 1] | —    | FFT of `Sig` (unshifted, complex-valued).                |
| `SigNewAdd`      | [M x 1] | —    | Newly arrived samples for the sliding-window update.     |
| `compset_rsb`    | [K x N] | —    | **Extracted component matrix**.                          |
| `Spec`           | [Nfft x 1] | — | Full single-sided amplitude spectrum from RSpectrendgene.|

### 1.4 VMD Parameters

| Variable     | Type    | Unit | Definition                                                                                |
|--------------|---------|------|-------------------------------------------------------------------------------------------|
| `alpha_vmd`  | scalar  | —    | **Bandwidth penalty parameter** (α). Larger values produce narrower-band modes.           |
| `tau_vmd`    | scalar  | —    | **Dual ascent time step**. 0 = exact reconstruction (noise-free). Default: **0**.         |
| `DC_vmd`     | scalar  | —    | Whether to fix the first mode at DC (0 Hz). 0 = no. Default: **0**.                       |
| `init_vmd`   | scalar  | —    | **Center frequency initialization**: 0 = all zeros, 1 = uniform spacing, 2 = random.      |
| `tol_vmd`    | scalar  | —    | **Convergence tolerance** for the ADMM loop. Default: **1×10⁻⁷**.                         |
| `K_vmd`      | scalar  | —    | **Number of modes** to decompose into. Default: **5**.                                    |
| `u_vmd`      | [K x N] | —    | Time-domain VMD modes (rows).                                                             |
| `kept_v`     | [k x N] | —    | Subset of modes whose Pearson correlation with `x_clean` exceeds `ce_rsb` (0.3).          |
| `time_vmd`   | scalar  | s    | Execution time of VMD.                                                                    |

### 1.5 EMD Parameters

| Variable     | Type    | Unit | Definition                                                                               |
|--------------|---------|------|------------------------------------------------------------------------------------------|
| `imfs`       | [N x K] | —    | **Intrinsic Mode Functions** from Empirical Mode Decomposition. Columns = modes.         |
| `nIMF`       | scalar  | —    | Number of extracted IMFs.                                                                |
| `kept_imf`   | [N x k] | —    | IMFs with Pearson correlation.                                                           |
| `time_emd`   | scalar  | s     | Execution time of EMD.                                                                  |

### 1.6 Bandpass Filtering Parameters

| Variable     | Type    | Unit | Definition                                                                               |
|--------------|---------|------|------------------------------------------------------------------------------------------|
| `sortInter`  | [B x 2] | Hz   | Sorted frequency band intervals `[fLow, fHigh]` from RSBAVMD spectrum trend analysis.     |
| `nBands`     | scalar  | —    | Number of frequency bands = `size(sortInter, 1)`.                                         |
| `fL`, `fH`   | scalar  | Hz   | Lower and upper cutoff frequencies for a bandpass sub-filter.                             |
| `bb`, `aa`   | vector  | —    | 4th-order Butterworth bandpass filter coefficients.                                       |
| `time_bp`    | scalar  | s     | Execution time of bandpass filtering.                                                     |
| `x_den_bp`   | [N x 1] | —    | **Bandpass denoised signal** = sum of filter outputs across all frequency bands.          |

### 1.7 Performance Metrics

| Variable       | Type           | Unit | Definition                                                                 |
|----------------|----------------|------|----------------------------------------------------------------------------|
| `m_MPE`        | scalar         | —    | **Embedding dimension** for Multiscale Permutation Entropy (MPE).          |
| `t_MPE`        | scalar         | —    | **Time delay** for phase space reconstruction.                             |
| `Scales`       | scalar         | —    | **Number of coarse-graining scales** for MPE.                              | 
| `snrComp(i)`   | scalar         | dB   | SNR of the i-th component relative to `x_clean`.                           |
| `corrComp(i)`  | scalar         | —    | Pearson correlation coefficient between the i-th component and `x_clean`.  |
| `mpeNoisy`     | scalar         | —    | MPE of the noisy signal `x_noisy`.                                         |
| `snrNoisy`     | scalar         | dB   | Original SNR = 2 dB (verification).                                        |
| `corrNoisy`    | scalar         | —    | Correlation between noisy signal and clean signal.                         |
| `mpeClean`     | scalar         | —    | MPE of the clean signal `x_clean` (reference).                             |

---

## 2. File: `demo2.m` — Measured Data RSBAVMD Decomposition

`demo2.m` applies RSBAVMD exclusively to **real measured vibration data** (from `data.mat`), visualizes the decomposed modes and their spectra, and compares the original and denoised signals.

---

## 3. Common Helper Functions (Both Files)

### 3.1 `RSBAVMD` — Recursive Sliding-window Band-Adaptive VMD

| Parameter    | Type         | Unit | Description                                                              |
|--------------|--------------|------|--------------------------------------------------------------------------|
| `Sig`        | [N x 1]      | —    | Time-domain input signal (column vector).                                |
| `SigFFT`     | [N x 1]      | —    | Previous frame FFT (complex, unshifted).                                 |
| `SigNewAdd`  | [M x 1]      | —    | New M samples for sliding window update.                                 |
| `SampFreq`   | scalar       | Hz   | Sampling frequency.                                                      |
| `tol`        | scalar       | —    | ADMM convergence tolerance.                                              |
| `ce`         | scalar       | —    | Correlation threshold (not used internally in RSBAVMD; passed through).  |
| `offset`     | scalar       | —    | Spectrum trend offset.                                                   |
| `cutpfreq`   | scalar       | —    | Normalized cutoff for Fourier fitting order.                             |
| `compset`    | [K x N]      | —    | Extracted component matrix.                                              |
| `Spec`       | vector       | —    | Amplitude spectrum from RSpectrendgene.                                  |

### 3.2 `RSpectrendgene` — Recursive Spectrum Trend Generation

| Output        | Type         | Unit | Description                                                                   |
|---------------|--------------|------|-------------------------------------------------------------------------------|
| `Spec`        | [N/2 x 1]   | —    | Single-sided amplitude spectrum (2·|FFT|/N).                                   |
| `Weight`      | [N/2 x 1]   | —    | Per-segment energy weight: mean squared spectrum within each frequency segment.|
| `SpecTrend`   | [N/2 x 1]   | —    | Smoothed spectrum envelope via Fourier basis fitting + offset.                 |
| `WeSpecTrend` | [N/2 x 1]   | —    | Weighted spectrum trend = `Weight .* SpecTrend`.                               |
| `sortInter`   | [B x 2]     | Hz   | Band intervals sorted by descending `WeightMax`. B = number of detected bands. |

### 3.3 `coef_ovefour` — Fourier Basis Fitting

Fits a linear combination of Fourier basis functions (cosine and sine) to the input spectrum for trend extraction.

| Parameter    | Type    | Definition                                                    |
|--------------|---------|---------------------------------------------------------------|
| `f`          | vector  | Input data (spectrum magnitudes).                             |
| `SampFreq`   | scalar  | Sampling frequency (for constructing the time axis).          |
| `orderAmp`   | scalar  | Half the number of Fourier basis functions.                   |
| `alpha`      | scalar  | Ridge regression regularization parameter. Fixed at **0.01**. |
| `fitF`       | vector  | Fitted values (the spectrum trend).                           |
| `finte`      | vector  | Integral of the fitted curve (used internally).               |

### 3.4 `RSVMD` — Recursive Sliding-window VMD

Implements VMD with incremental FFT for sliding-window processing. Solves the variational problem via ADMM.

| Parameter    | Type         | Unit | Description                                                                        |
|--------------|--------------|------|------------------------------------------------------------------------------------|
| `Alpha`      | scalar/[1,K] | —    | Bandwidth penalty factor(s). One per mode, or scalar broadcast to all modes.       |
| `tau`        | scalar       | —    | Dual ascent step size. 0 = exact reconstruction.                                   |
| `K`          | scalar       | —    | Number of modes.                                                                   |
| `DC`         | scalar       | —    | Fix first mode at DC (0=no, 1=yes).                                                |
| `init`       | scalar       | —    | Center frequency initialization: 0=zeros, 1=uniform, 2=random.                     |
| `tol`        | scalar       | —    | Convergence tolerance.                                                             |
| `M`          | scalar       | samples | Sliding step size = length(SigNewAdd). 0 = cold start (no sliding).             |
| `u`          | [K x N]      | —    | Time-domain modes.                                                                 |
| `u_hat`      | [N x K]      | —    | Mode spectra (complex, unshifted).                                                 |
| `lambda_hat` | [N x 1]      | —    | Lagrangian multiplier (dual variable) in frequency domain.                         |

### 3.5 `recursive_FFT` — Incremental FFT Update

| Parameter    | Type       | Unit | Definition                                                                   |
|--------------|------------|------|-----------------------------------------------------------------------------|
| `sigOrigin`  | [Nf x 1]   | —    | Previous frame time-domain signal.                                           |
| `sigFFT`     | [Nf x 1]   | —    | Previous frame FFT (complex, unshifted).                                     |
| `sigNewAdd`  | [k x 1]    | —    | New k samples replacing the oldest k samples.                                |

### 3.6 `updateCenterFreq` — Spectral Centroid

Computes the weighted spectral centroid (center of gravity) of a mode's spectrum, used to update `omega_k` in the VMD ADMM loop.

### 3.7 `VMD` — Standard Variational Mode Decomposition

Classical VMD implementation with signal mirroring extension. Identical parameter semantics to `RSVMD` but without sliding-window support.

### 3.8 `MPerm` — Multiscale Permutation Entropy (MPE)

| Parameter    | Type    | Definition                                                     |
|--------------|---------|----------------------------------------------------------------|
| `X`          | vector  | Input time series.                                             |
| `MPE`        | scalar  | Normalized MPE value. Normalized by `log₂(m!)`. Range: [0, 1]. |

## 4. External Dependencies

| File / Toolbox          | Purpose                                                       |
|-------------------------|---------------------------------------------------------------|
| `data.mat`              | Contains real measured vibration data (`Sig` field).           |
| MATLAB Signal Processing Toolbox | `butter`, `filtfilt`, `hilbert` functions.            |
| MATLAB `emd` function   | Empirical Mode Decomposition (available in Signal Processing Toolbox or third-party packages). |

---

## 5. Contextual Notes

1. **Sliding window mechanism**: Both scripts simulate online/streaming processing by using a window step size of M=50 samples.

2. **Frequency band detection**: `RSpectrendgene` uses negative peak detection on the spectrum trend curve to partition the frequency axis. 

3. **Alpha adaptation**: In RSBAVMD, the bandwidth penalty `α_i` is scaled inversely with each band's relative bandwidth: narrow bands get larger α (tighter frequency constraint), wide bands get smaller α (more flexibility). The base value β=500 was empirically chosen.
