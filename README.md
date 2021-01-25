# FIR_to_IIR
Supplementary for resamplers

Attempt to convert long FIR kernel for resamplers to smaller length IIR kernel. The main task is to find such coefficients if IIR filter that impulse response match as good as possible the FIR filter kernel. But the number of IIR filter coefficients must be smaller in compare of total number of FIR filter coefficients.
