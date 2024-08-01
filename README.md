# Skew-Spectra
Optimal Extraction of the Non-Gaussian Information of the Cosmic Large-scale Structure

Multiple scripts written during my 2024 internship on "Optimal Extraction of the Non-Gaussian Information of the Cosmic Large-scale Structure" at the LAPth at Annecy supervised by Azadeh Moradinezhad Dizgah.

Mathematica files are responsible for generating expressions necessary for calculating the expectation value of the Bi-spectrum at a given momentum k.

randomrange files are C code that aim to find an optimal number of biases and their values, that minimise the number of UV/IR divergences needed after fftlog calculations

skew_spectra_fftlog C code calculates the expectation value of the Bi-spectrum for a range of k momenta values.


Known Issues:

- Line of sight integral functions take time to calculate, especially higher orders. So far orders 1 to 10 have been calculated, orders 11 and 12 are necessary for some terms and still need to be generated.

- The UV/IR divergences require limits that are slow to perform. Either optimising the current implementation or finding a faster method is necessary.

Sources:

Bispectrum decomposition: https://arxiv.org/pdf/2010.14267

fftlog method: https://arxiv.org/abs/1708.08130