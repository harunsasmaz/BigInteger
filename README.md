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

* You will see [Garner's Algorithm](https://cp-algorithms.com/algebra/chinese-remainder-theorem.html#toc-tgt-2) under <code>gcd.h</code>, it is a corollary of Chinese Remainder Theorem.

* Under <code>gcd.h</code>, there are two different GCD algorithms, [euclidian](https://en.wikipedia.org/wiki/Euclidean_algorithm) and [binary](https://en.wikipedia.org/wiki/Binary_GCD_algorithm). You are expected to understand these GCD algorithms and their inverse functions.

* Primality tests of long precision numbers are pretty hard with brute-force way. Hence, this implementation uses [Miller-Rabin Algorithm](https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test) to check if a number is prime or not.


## Test

* You are given a pre-defined test class which simply tests operations of BigInteger class. You may add other test criterias if you wish or you can simply run <code>main.cpp</code> to see results.
