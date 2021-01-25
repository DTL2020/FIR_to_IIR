# FIR_to_IIR
Supplementary for resamplers

Attempt to convert long FIR kernel for resamplers to smaller length IIR kernel. The main task is to find such coefficients if IIR filter that impulse response match as good as possible the FIR filter kernel. But the number of IIR filter coefficients must be smaller in compare of total number of FIR filter coefficients.

Current version of solver works about good at iM from about 3 to 6 but for unknown reason starts to fails if iM > 6..7. May be because of rounding errors with even double calculations in linear system solver.
