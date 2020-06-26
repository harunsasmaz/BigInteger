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

BN::BN(const string &str, RadixBI radix)
{
    if(radix == RadixBI::dec)
    {
        BN bn(1,0);
        for (auto i : str) {
            if (i < '0' || i > '9')
                continue;
            bn.mulbaseappr(10);
            bn += BN(i - '0');
        }
        data = move(bn.data);
        return;
    }

    // radix: hex
    size_t length = str.size();
    data.resize((length + bz * 2 - 1) / (bz * 2));

    size_t index = data.size() - 1;
    size_t shift = bz * 2 - (data.size() * bz * 2 - length);

    for (size_t i = 0; i < length; ++i) {
        bt d;
        if (str[i] >= '0' && str[i] <= '9')
            d = str[i] - '0';
        else if (str[i] >= 'A' && str[i] <= 'F')
            d = str[i] - 'A' + 10;
        else if (str[i] >= 'a' && str[i] <= 'f')
            d = str[i] - 'a' + 10;
        else
            throw invalid_argument(string("BN constructor: invalid char '") + str[i] + "' in HEX string");
        data[index] = (data[index] << 4) | d;
        --shift;
        if (shift == 0) {
            shift = bz * 2;
            --index;
        }
    }
    norm(data);
}

// Swap operation

void BN::swap(BN &bn) noexcept
{
    data.swap(bn.data);
}

// Assignment operations

BN & BN::operator = (const BN& bn)
{
    if (this == &bn)
        return *this;

    data = bn.data;
    return *this;
}

BN & BN::operator = (BN&& bn) noexcept
{
    if (this == &bn)
        return *this;

    data = move(bn.data);
    return *this;
}

// Addition operations

const BN BN::operator + (const BN& bn) const 
{
    return (move(BN(*this)) += bn);
}

BN & BN::operator ++()
{
    for (auto& i : data) {
        ++i;
        if (i != 0)
            return *this;
    }
    data.emplace_back(1);
    return *this;
}

BN& BN::operator += (const BN& bn) {
    data.resize(max(data.size(), bn.data.size()));

    bt2 sum = 0;
    for(size_t pos = 0; pos < bn.data.size(); pos++) {
        sum = (sum >> bz8) + data[pos] + bn.data[pos];
        data[pos] = sum;
    }

    for(size_t pos = bn.data.size(); pos < data.size(); pos++) {
        sum = (sum >> bz8) + data[pos];
        data[pos] = sum;
    }

    sum >>= bz8;
    if (sum)
        data.emplace_back(sum);
    return *this;
}

// Subtraction operations

const BN BN::operator - (const BN& bn) const
{
    return move(BN(*this) -= bn);
}

BN & BN::operator -- () noexcept
{
    for (auto& i : data) {
        if (i) {
            --i;
            return *this;
        }
        --i;
    }

    norm(data);
    return *this;
}

BN& BN::operator -= (const BN& bn)
{
    if (data.size() < bn.data.size())
        data.resize(bn.data.size());

    bool flag = 0;
    size_t pos = 0;
    for (; pos < bn.data.size() && pos < data.size(); ++pos) {
        bt2s res = static_cast<bt2s>(data[pos]) - bn.data[pos] - flag;
        data[pos] = static_cast<bt>(res);
        flag = (res < 0);
    }

    for (; flag && pos < data.size(); ++pos) {
        if (data[pos])
            flag = false;
        --data[pos];
    }

    norm(data);
    return *this;
}

// Base operations

BN BN::mulbt(size_t t) const
{
    if(t == 0)
        return *this;

    BN res(data.size() + t, 0);
    for(size_t i = 0; i < data.size(); ++i)
        res.data[i + t] = data[i];
    norm(res.data);
    return res;
}

BN BN::divbt(size_t t) const
{
    if(t >= data.size())
        return BN::BigInteger0();

    return move(vector<bt>(data.begin() + t, data.end()));
}

BN BN::modbt(size_t t) const
{
    if(t >= data.size())
        return *this;

    return move(vector<bt>(data.begin(), data.begin() + t));
}

const BN BN::mulbase(const bt & multiplier) const
{
    return move(BN(*this).mulbaseappr(multiplier));
}

BN& BN::mulbaseappr(const bt &multiplier)
{
    if (!multiplier) {
        data.resize(1);
        data.front() = 0;
        return *this;
    }

    bt2 curr = 0;
    for(auto& i : data) {
        i = curr += static_cast<bt2>(i) * multiplier;
        curr >>= bz8;
    }
    if (curr)
        data.emplace_back(curr);
    return *this;
}

const BN BN::divbase(const bt & divider) const 
{
    return move(BN(*this).divbaseappr(divider));
}

BN& BN::divbaseappr(const bt &diviser)
{
    if(diviser == 0)
        throw divide_by_zero();

    bt2 curr = 0;
    for (size_t i = data.size(); i; --i) {
        curr = (curr % diviser << bz8) + data[i - 1];
        data[i - 1] = curr / diviser;
    }

    norm(data);
    return *this;
}

const BN BN::modbase(const bt &diviser)const
{
    if(diviser == 0)
        throw divide_by_zero();

    bt2 curr=0;
    for (size_t i = data.size()- 1; i < data.size(); --i) {
        curr <<= bz8;
        curr += data[i];
        curr %= diviser;
    }

    BN result;
    result.data[0]=curr;
    return result;
}

BN& BN::modbaseappr(const bt &diviser)
{
    if(diviser == 0)
        throw divide_by_zero();
    bt2 curr = 0;
    for(size_t i = data.size() - 1; i < data.size(); --i) {
        curr <<= bz8;
        curr += data[i];
        data[i] = curr / diviser;
        curr %= diviser;
    }
    
    data.resize(1);
    data[0] = curr;
    return *this;
}