# ELEC9123 Design Task B - BPSK Wireless Communication System Analysis

## Project Information
- **Course**: ELEC9123 Digital Communications
- **Task**: Design Task B - BPSK Performance Analysis
- **Author**: [Your Name]
- **Student ID**: [Your zID]
- **Date**: 2025

## Overview
This project implements comprehensive BPSK (Binary Phase Shift Keying) communication system performance analysis across multiple levels, from basic noise channel analysis to advanced fading scenarios with random deployment.

## System Requirements
- **MATLAB Version**: R2023b or later
- **Required Toolboxes**: Communications Toolbox, Statistics and Machine Learning Toolbox
- **Memory**: 8GB RAM recommended for Level 5
- **Processing**: Multi-core processor recommended

## File Structure
```
level2.m                 # Level 2: AWGN + Laplacian noise analysis
level5.m                 # Level 5: Advanced fading channels + random deployment
README.md                # This documentation file
```

## Implementation Levels

### Level 2 (Pass Level) - `level2.m`
**Channel Models:**
- AWGN (Additive White Gaussian Noise)
- Laplacian impulsive noise

**Key Features:**
- Custom Laplacian random number generator using inverse transform method
- BER simulation with 10⁷ bits per SNR point
- Theoretical validation: AWGN P_b = Q(√(2E_b/N₀)), Laplacian P_b = (1/2)exp(-√(2E_b/N₀))
- Statistical verification of Laplacian distribution (PDF/CDF)

**Outputs:**
- BER performance comparison plot
- Laplacian distribution verification plots

### Level 5 (High Distinction) - `level5.m`
**Channel Models:**
1. AWGN baseline
2. Laplacian noise
3. Rayleigh fading + AWGN (fixed distance d=1m)
4. Rayleigh fading + AWGN + random deployment
5. Rician fading + AWGN (K-factor = 5, fixed distance)
6. Rician fading + AWGN + random deployment

**Advanced Features:**
- Monte Carlo simulation with 10⁶ samples per SNR point
- Random deployment modeling with path-loss exponent α = 2.2
- Outage probability analysis with rate threshold C = 1.2 bps/Hz
- Numerically stable Marcum Q-function implementation
- Advanced integration techniques for analytical calculations

**Outputs:**
- BER performance in noise-only channels
- BER performance in fading channels
- Outage probability analysis

## How to Run

### Level 2 Execution
```matlab
level2
```
**Runtime:** ~2-5 minutes
**Memory:** ~1GB

### Level 5 Execution
```matlab
zID_LastName_DTB_2025_Corrected
```
**Runtime:** ~15-30 minutes
**Memory:** ~2-4GB

## Key Technical Implementations

### Laplacian Random Number Generation
```matlab
function r = laprnd(m,n,mu,b)
    U = rand(m,n) - 0.5;
    r = mu - b * sign(U) .* log(1 - 2*abs(U));
end
```

### Stable Marcum Q-Function (Level 5)
```matlab
function q = MarcumQ1_stable(a, b)
    integrand = @(t) t .* besseli(0, a.*t, 1) .* exp(-(t-a).^2 / 2);
    q = integral(integrand, b, inf);
end
```

### Random Deployment Modeling (Level 5)
```matlab
dRand = R*sqrt(rand(N,1));           % Uniform in circular area
path_loss = dRand.^alpha;            % Path-loss with exponent α
effective_snr = g_lin * |h|² / path_loss;
```

## Performance Metrics

### BER Analysis
- **Simulation Method**: Monte Carlo with hard decision detection
- **Theoretical Formulas**: Verified against established communication theory
- **SNR Range**: 0 to 15 dB (E_b/N₀)

### Outage Probability (Level 5)
- **Threshold**: γ_th = 2^C - 1 where C = 1.2 bps/Hz
- **Analysis**: Both fixed distance and random deployment scenarios
- **Integration**: Analytical calculation with numerical verification

## Expected Results

### Channel Performance Ranking (Best to Worst BER)
1. AWGN
2. Rician fading (fixed distance)
3. Laplacian noise
4. Rayleigh fading (fixed distance)
5. Rician fading (random deployment)
6. Rayleigh fading (random deployment)

### Key Observations
- Random deployment significantly degrades performance due to path-loss
- Rician channels outperform Rayleigh due to line-of-sight component
- Laplacian noise shows exponential BER decay
- Outage probability increases with random deployment

## Validation and Accuracy

### Statistical Validation
- **Level 2**: Large sample size (10⁷ bits) ensures high accuracy
- **Level 5**: Monte Carlo error O(1/√N) with N = 10⁶
- **Theory Matching**: Simulation results validate analytical formulas

### Numerical Stability
- Robust implementations prevent overflow/underflow
- Stable integration algorithms for complex scenarios
- Error handling for edge cases

## Troubleshooting

### Common Issues
1. **Memory Limitations**: Reduce sample size or process sequentially
2. **Long Runtime**: Consider parallel processing
3. **Integration Convergence**: Adjust tolerance parameters
4. **Figure Display**: Ensure graphics drivers are updated

### Performance Tips
- Monitor memory usage during Level 5 execution
- Use vectorized operations for efficiency
- Verify numerical stability at extreme SNR values

## Requirements Compliance

### Level 2 ✅
- [x] AWGN channel BER analysis
- [x] Laplacian noise implementation
- [x] Custom random number generator
- [x] Statistical distribution verification
- [x] Performance comparison plots

### Level 5 ✅
- [x] All Level 2 requirements
- [x] Rayleigh and Rician fading channels
- [x] Random deployment scenarios
- [x] Outage probability analysis
- [x] Advanced numerical implementations
- [x] Comprehensive performance validation

## References
1. Proakis, J. G., & Salehi, M. (2008). Digital Communications (5th ed.)
2. Simon, M. K., & Alouini, M. S. (2005). Digital Communication over Fading Channels
3. Goldsmith, A. (2005). Wireless Communications
4. MATLAB Documentation: Communications Toolbox

## Contact Information
- **Student**: [Your Name] ([Your zID]@student.unsw.edu.au)
- **Course**: ELEC9123 Digital Communications
- **Institution**: UNSW Sydney

---
*Complete implementation for ELEC9123 Design Task B - Levels 2 and 5*
