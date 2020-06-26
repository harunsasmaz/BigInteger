#include "BigInteger.h"

#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stack>
#include <stdexcept>

using namespace std;

#define BN BigInteger

namespace {

    constexpr bt4 max_size_fast_mul = numeric_limits<bt4>::max() / bmax / bmax / bmax * (bmax - 1) - 1;

    constexpr size_t karatsuba_min_size = 50;

    inline void norm(vector<bt> & ba) noexcept
    {
        while(ba.size() > 1 && ba.back() == 0)
            ba.pop_back();
    }

    const BN reductionBarrettPrecomputation(const BN& mod) {
        size_t rbc = mod.digit_count() * 2 + 1;
        vector<bt> ba(rbc);
        ba.back() = 1;
        return move(BN(ba, rbc) / mod);
    }

    constexpr size_t KarySize = 256;
    constexpr bt KaryMask = 0xFF;
    constexpr size_t KaryBits = 8;

} // namespace

// Divide by zero logical arithmetic error struct

divide_by_zero::divide_by_zero() : logic_error("Divide by Zero") {}

// CONSTRUCTORS 

BN::BN() : data(1) {}

BN::BN(uint64_t base, int) : data(base) {} 

BN::BN(const BN& bn) : data(bn.data) {}

BN::BN(BN&& bn) noexcept : data(move(bn.data)) {}

BN::BN(vector<bt> && v) noexcept : data(move(v)) { norm(data); }

BN::BN(const vector<bt> & v) : data(v) { norm(data); }

BN::BN(vector<bt> && v, size_t rbc) : data(move(v)) { data.resize(rbc); }

BN::BN(const vector<bt> & v, size_t rbc) : data(v.begin(), v.begin() + rbc) {}

BN::BN(uint64_t x) : data((sizeof(uint64_t) + bz - 1) / bz)
{
    for(size_t i = 0; i < data.size(); i++, x >>= bz8)
        data[i] = static_cast<bt>(x);
    norm(data);
}


