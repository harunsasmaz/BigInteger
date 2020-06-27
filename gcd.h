#ifndef _GCD_H
#define _GCD_H

#include "BigInteger.h"
#include "BigIntegerSign.h"

#include <vector>

BN gcdEuclidean(BN x,BN y);
BN gcdBinary(BN x,BN y);

BN gcdInverseEuclidean(const BN& a, const BN& mod);
BN gcdInverseEuclideanBinary(BN x,BN mod);

std::vector <BN> multi_inverse(const std::vector <BN>& x, const BN &mod);
BN Garner(const std::vector <BN>& m, const std::vector <BN>& v);
BN CTO(const std::vector <BN>& m, const std::vector <BN>& v);

bool isPrimeMillerRabin(const BN& n, size_t k);

#endif