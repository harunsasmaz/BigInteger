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

// Base multiplciation
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

// Base division
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

// Base modulo
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

// Multiplication operations

const BN BN::operator * (const BN& bn) const
{
    return karatsuba_mul(bn);
}

const BN BN::naive_mul(const BN& bn) const {
    // classical O(n*n) multiplication.
    // b * a is slightly faster than a * b
    const BN& b = data.size() > bn.data.size() ? *this : bn;
    const BN& a = data.size() > bn.data.size() ? bn : *this;

    if (a.data.size() == 1)
        return b.mulbase(a.data.front());


    BN result(a.data.size() + b.data.size(), 0);
    for (size_t i = 0; i < b.data.size(); ++i) {
        bt2 curr = 0;
        bt2 x = b.data[i];
        for (size_t j = 0; j < a.data.size(); ++j) {
            curr = (curr >> bz8) + result.data[i + j] + x * a.data[j];
            result.data[i + j] = curr;
        }
        result.data[i + a.data.size()] = curr >> bz8;
    }
    norm(result.data);
    return result;
}

const BN BN::fast_mul(const BN& bn) const {
    size_t n = data.size();
    size_t m = bn.data.size();

    BN result(n + m, 0);

    bt4 t = 0;
    for(size_t s = 0; s < m + n - 1; s++) {

        size_t end_index = min(n - 1, s);
        size_t start_index = s >= m ? s - m + 1 : 0;
        for(size_t i = start_index, j = s - start_index; i <= end_index; i++, --j)
            t += static_cast<bt2>(data[i]) * bn.data[j];


        result.data[s] = t;
        t = t >> bz8;
    }

    result.data.back() = t;
    norm(result.data);
    return move(result);
}

/* TODO: Karatsuba Multiplication Implementation */

/* TODO: divmod function implementation */

const BN BN::operator / (const BN& bn) const
{
    BN div,mod;
    divmod(bn,div,mod);
    return move(div);
}

const BN BN::operator % (const BN& bn) const
{
    BN div,mod;
    divmod(bn,div,mod);
    return move(mod);   
}

// Conditional operations

bool BN::operator < (const BN& bn) const noexcept
{
    // TODO: maybe can be replaced to compare ba and bn.ba
    if(data.size() > bn.data.size())
        return false;
    if(data.size() < bn.data.size())
        return true;
    for(size_t i = data.size() - 1; i < data.size(); i--)
    {
        if(data[i] > bn.data[i])
            return false;
        if(data[i] < bn.data[i])
            return true;
    }

    return false;
}

bool BN::operator > (const BN& bn) const noexcept
{
    return bn < *this;
}

bool BN::operator <= (const BN& bn) const noexcept
{
    return !(*this > bn);
}

bool BN::operator >= (const BN& bn) const noexcept
{
    return !(*this < bn);
}

bool BN::operator == (const BN& bn) const noexcept
{
    return data == bn.data;
}

bool BN::operator != (const BN& bn) const noexcept
{
    return !(*this == bn);
}

bt BN::operator [](size_t index) const noexcept
{
    return data[index];
}

// Digit information

size_t BN::digit_count() const noexcept
{
    return data.size();
}

size_t BN::bit_count() const noexcept
{
    size_t x = 0;
    bt value = data.back();
    while (value) {
        ++x;
        value >>= 1;
    }
    return (data.size() - 1) * bz8 + x;
}

// Shift operations

const BN BN::operator >> (size_t shift) const {
    return move(BN(*this) >>= shift);
}

const BN BN::operator << (size_t shift) const {
    return move(BN(*this) <<= shift);
}