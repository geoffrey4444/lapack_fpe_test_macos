# lapack_fpe_test_macos

This code tests for floating-point exceptions (overflow, invalid operation, divide by zero) 
when solving a generalized eigenvalue problem on macOS. 

There are at least two ways to compile this code. The first uses Apple's Accelerate framework:

<code>clang++ -framework Accelerate -o Test_dggev_Accelerate Test_dggev.cpp</code>

A second way uses homebrew (https://brew.sh) to install openblas:

<code>brew install openblas</code>

<code>clang++ -L/usr/local/opt/openblas/lib -lopenblas -o Test_dggev Test_dggev.cpp</code>

There should be no floating-point exception, in which case the code prints the real part of the first eigenvalue found. 
If there is a floating-point exception, the code will terminate before printing any results.

If you have any questions about this repo, please contact Geoffrey Lovelace 
(glovelace at fullerton dot edu).
