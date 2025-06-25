# BP-decoder-for-NB_LDPC-codes Matlab with C mex files

This repository contains a platform for performing Non-Binary Low-Density Parity-Check (NB-LDPC) sum-product decoding, implementing both flooding and layered approaches. It leverages the Fast Fourier Transform-QSPA (FFT-QSPA) algorithm developed by MacKay and Davey (Chapter 3.4 [1]), utilizing 2,3 values in the first parameter. Additionally, it includes majority decoding with 3,4 values. In commented code done Binary image optimization using weight spectrum.



*   **NB-LDPC Decoding:** Implements sum-product decoding for Non-Binary LDPC codes.
*   **Flooding and Layered Decoding:** Supports both flooding and layered decoding algorithms.
*   **FFT-QSPA Implementation:** Utilizes the FFT-QSPA algorithm developed by MacKay and Davey, enabling efficient decoding.
*   **Prime Size Flexibility:** The use of FFT allows for different prime sizes up to 1024.
*   **Majority Decoding:** Includes majority decoding with 3,4 values.

**References:**

1.  [M. C. Davey, *Error-correction using low-density parity-check codes*, Ph.D. dissertation, Univ. of Cambridge, Cambridge, U.K., Dec. 1999.](https://github.com/Lcrypto/BP-decoder-for-NB_LDPC-codes/blob/master/davey_phd.pdf) 
2.  D. Declercq and M. Fossorier, "Decoding Algorithms for Nonbinary LDPC Codes Over GF (q) ," in *IEEE Transactions on Communications*, vol. 55, no. 4, pp. 633-643, April 2007.

