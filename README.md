# BigInteger

A C++ implementation of BigInteger a.k.a BigNumber class

## Before You Start

* BigNumber class has a Cryptographic usage and also used in Crypto Currencies.

    - For Cryptography [check](https://case.edu/affil/sigmaxi/files/CryptoslidesSinger.pdf)
    - For Crypto Currencies (Ethereum) [check](https://docs.ethers.io/v5/api/utils/bignumber/)

* I strongly suggest you to learn [string](http://www.cplusplus.com/reference/string/string/)  and [vector](http://www.cplusplus.com/reference/vector/vector/) modules of C++

* Rotation operations use [K-ary Trees](https://xlinux.nist.gov/dads/HTML/karyTree.html). Please, understand how left-right rotations are made in this data structure

* This class makes use of [Karatsuba Algorithm](https://en.wikipedia.org/wiki/Karatsuba_algorithm#Pseudocode)
So, it is better for you to have a look at how this algorithm works

* This class makes use of [Barret's Reduction Algorithm](https://en.wikipedia.org/wiki/Barrett_reduction). So, please also have a look at how this algorithm works.
