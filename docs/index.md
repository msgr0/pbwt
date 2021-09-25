# journey in the wild(card)pbwt
---
## bitvector
The bitvector is used to calculate in fixed time whether or not an open blocks ends at the current `k` index.
The idea of the bitvector for wildcards is to build, for each allele, a vector of 0 and 1 meaning a match in k+1 for a\[i\] and his predecessor a\[i-1\].
We have t-bitvectors in order to anticipately correct a wildcard with the respectiv t-allele.

