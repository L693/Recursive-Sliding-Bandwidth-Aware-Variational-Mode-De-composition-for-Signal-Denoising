
1. Dependencies

1.1 MATLAB Version Requirements

Minimum Compatible Version: MATLAB R2019a

Recommended Version: MATLAB R2021a and above (faster execution speed, better figure display quality)

1.2 Required Toolboxes

Signal Processing Toolbox and Wavelet Toolbox

1.3 No Additional Third-Party Dependencies

All core algorithms (RSBAVMD, RSVMD, recursive FFT, Weighted Spectrum Trend) are implemented in pure MATLAB

EEMD uses a built-in compatible implementation (eemd_compatible.m), no need to install additional EEMD toolboxes

All evaluation metric functions are custom-implemented and do not depend on any external libraries

2. Project Structure

2.1 Recommended Standard Repository Structure

RSBAVMD-Signal-Denoising/

├── main.m                          # Main experiment script (your complete code)

├── core/                           # Core algorithm modules

│   ├── RSBAVMD.m                   # RSBAVMD main algorithm entry

│   ├── RSVMD.m                     # Recursive VMD core algorithm

│   ├── VMD.m                       # Standard VMD algorithm (consistent with the paper)

│   └── RSpectrendgene.m            # Weighted Spectrum Trend (WST) spectral partitioning module

├── utils/                          # General utility functions

│   ├── calculate_metrics.m         # All evaluation metric calculation functions

│   ├── eemd_compatible.m           # Cross-version compatible EEMD implementation

│   ├── recursive_FFT.m             # Recursive sliding Fourier transform

│   ├── fft_spectrum.m              # Standardized FFT spectrum calculation

│   └── coef_ovefour.m              # Overcomplete Fourier fitting function

├── config.m                        # Global parameter configuration file

├── data/                           # Data directory

│   ├── synthetic/                  # Synthetic signal generation scripts

│   └── real/                       # Real vortex-induced vibration signal storage directory

├── results/                        # Result output directory

│   ├── figures/                    # Generated figure files

│   └── metrics/                    # Quantitative metric data files

└── README.md                       # Project documentation

2.2 Module Division of Your Current Single-File Code

The code already contains all the above functions, logically divided into 7 modules:

1. Synthetic Signal Generation Module: Generates composite test signals containing damped vibration and amplitude modulation components.

2. Evaluation Metrics Module: SNR, MSE, RMSE, correlation coefficient, envelope smoothness, permutation entropy, multiscale permutation entropy.

3. Comparison Algorithms Module: EMD, EEMD, bandpass filtering, grid-optimized VMD, wavelet denoising.

4. RSBAVMD Core Module: Recursive spectrum trend generation, adaptive parameter selection, recursive VMD decomposition.

5. Parameter Sensitivity Analysis Module: Effect of window length and sliding step on algorithm performance.

6. SNR Robustness Analysis Module: Algorithm performance comparison under different input signal-to-noise ratios.

7. Real Signal Validation Module: Denoising effect comparison of measured vortex-induced vibration signals.
