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

// Karatsuba Multiplication Algorithm

inline vector<bt> karatsubaSum(const vector<bt>& u, size_t start, size_t n, size_t m) 
{
    vector<bt> result;
    result.reserve(m + 1);

    bt2 sum = 0;
    for (size_t pos = 0; pos < n; ++pos) {
        sum += u[start + pos] + u[start + pos + n];
        result.emplace_back(sum);
        sum >>= bz8;
    }
    if (n != m) {
        sum += u[start + n + n];
        result.emplace_back(sum);
        sum >>= bz8;
    }
    result.emplace_back(sum);
    return move(result);
}

vector<bt> karatsubaRecursive(const vector<bt>& U, const vector<bt>& V, size_t start, size_t count) 
{
    const size_t len = count, n = len / 2, m = len - n;

    if (n < karatsuba_min_size) {
        vector<bt> result;
        result.reserve(len + len);

        bt4 t = 0;
        for(size_t s = 0; s < len + len - 1; s++) {

            size_t end_index = min(len - 1, s) + start;
            size_t start_index = s >= len ? s - len + 1 : 0;
            for(size_t i = start_index + start, j = s - start_index + start; i <= end_index; i++, --j)
                t += static_cast<bt2>(U[i]) * V[j];


            result.emplace_back(t);
            t = t >> bz8;
        }
        result.emplace_back(t);
        return move(result);
    }

    const vector<bt>& u01 = karatsubaSum(U, start, n, m);
    const vector<bt>& v01 = karatsubaSum(V, start, n, m);

    vector<bt> A = move(karatsubaRecursive(U, V, start + n, m));
    vector<bt> B = move(karatsubaRecursive(U, V, start, n));
    const vector<bt> C = move(karatsubaRecursive(u01, v01, 0, m + 1));

    vector<bt> result = B;
    result.resize(len + len);

    for (size_t i = 0; i < A.size(); ++i)
        result[i + n + n] = A[i];

    const size_t abcSize = m + m + 2;
    A.resize(abcSize);
    B.resize(abcSize);

    bt2s sum = 0;
    for (size_t i = 0; i < abcSize; ++i) {
        sum += result[i + n];
        sum += C[i];
        sum -= A[i];
        sum -= B[i];
        result[i + n] = sum;
        sum >>= bz8;
    }
    for (size_t i = n + abcSize; i < result.size(); ++i) {
        sum += result[i];
        result[i] = sum;
        sum >>= bz8;
    }
    return move(result);
}

const BN BN::karatsuba_mul(const BN& bn) const 
{
    size_t x = data.size(), y = bn.data.size(), len = max(x, y);

    if(min(x, y) < karatsuba_min_size)
        return this->fast_mul(bn);

    vector<bt> U = data;
    vector<bt> V = bn.data;
    U.resize(len); V.resize(len);

    return karatsubaRecursive(U, V, 0, len);
}

// divmod operation and division operators

void BN::divmod(const BN& bn, BN& div, BN& mod) const
{
    if(bn.isZero())
        throw divide_by_zero();

    if(bn.data.size() == 1) {
        div = move(this -> divbase(bn.data.front()));
        mod = move(this -> modbase(bn.data.front()));
        return;
    }

    if(*this < bn) {
        div = BN::BigInteger0();
        mod = *this;
        return;
    }

    BN dividend(*this), divider(bn);

    bt d = bsize  / (bn.data.back() + 1);
    if (d != 1) {
        dividend.mulbaseappr(d);
        divider.mulbaseappr(d);
    }

    size_t n = divider.data.size();
    size_t m = dividend.data.size() - n + 1;

    dividend.data.resize(dividend.data.size() + 2);
    divider.data.resize(divider.data.size() + 1);

    div.data.resize(m + 1);

    vector<bt> temp(n + 1);
    for (size_t j = m; j <= m; --j) {
        bt2 q = (dividend.data[j + n] * bsize + dividend.data[j + n - 1]) / divider.data[n-1];
        bt2 r = (dividend.data[j + n] * bsize + dividend.data[j + n - 1]) % divider.data[n-1];

        if (q == bsize || q * divider.data[n-2] > bsize * r + dividend.data[j + n - 2]) {
            --q;
            r += divider.data[n-1];
            if (r < bsize && q * divider.data[n-2] > bsize * r + dividend.data[j + n - 2])
                --q;
        }

        if (!q) {
            div.data[j] = 0;
            continue;
        }

        bt4s x = 0;
        for (size_t i = 0; i < n; ++i) {
            x += dividend.data[j + i];
            x -= q * divider.data[i];
            dividend.data[j + i] = x;
            x >>= bz8;
        }
        x += dividend.data[j + n];
        dividend.data[j + n] = x;

        // If `x' is negative, than `q' is too large.
        // Decrement `q' and update `dividend'.
        if (x < 0) {
            --q;
            x = 0;
            for (size_t i = 0; i < n; ++i) {
                x += dividend.data[j + i];
                x += divider.data[i];
                dividend.data[j + i] = x;
                x >>= bz8;
            }
            x += dividend.data[j + n];
            dividend.data[j + n] = x;
        }

        div.data[j] = q;
    }

    norm(div.data);
    norm(dividend.data);

    if (d != 1)
        dividend.divbaseappr(d);
    mod = move(dividend);
}

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

BN& BN::operator >>= (size_t shift) {
    if (shift == 0)
        return *this;
    size_t baseshift = shift / bz8;
    size_t realshift = shift - baseshift * bz8;

    if (realshift == 0)
        return *this = divbt(baseshift);

    if (baseshift >= data.size()) {
        data.resize(0);
        data[0] = 0;
        return *this;
    }

    for (size_t i = 0; i < data.size() - baseshift - 1; ++i) {
        data[i] =
            (data[i + baseshift] >> realshift) |
            (data[i + baseshift + 1] << (bz8 - realshift));
    }
    data[data.size() - baseshift - 1] = data.back() >> realshift;
    data.resize(data.size() - baseshift);

    norm(data);
    return *this;
}

BN& BN::operator <<= (size_t shift) {
    if(shift == 0)
        return *this;

    size_t baseshift = shift / bz8;
    size_t realshift = shift - baseshift * bz8;

    if (realshift == 0)
        return *this = mulbt(baseshift);

    data.resize(data.size() + baseshift + 1, 0);
    data.back() = data.back() >> (bz8 - realshift);
    for (size_t i = data.size() - 1; i; --i) {
        data[i + baseshift] =
            (data[i - 1] >> (bz8 - realshift)) |
            (data[i] << realshift);
    }
    data[baseshift] = data[0] << realshift;
    norm(data);
    return *this;
}

// Power operations

BN BN::reductionBarrett(const BN& mod, const BN& mu) const {

    size_t k = mod.data.size();
    if(k * 2 < data.size()) {
        return (*this) % mod;
    }

    BN q1 = divbt(k-1), q2 = q1*mu, q3 = q2.divbt(k+1);
    BN r1 = modbt(k + 1), r2 = (q3 * mod).modbt(k+1);
    r1 -= r2;
    while (r1 >= mod)
        r1 -= mod;
    return r1;
}

BN BN::pow(uint64_t power) const
{
    if (power == 0)
        return BN::BigInteger1();

    BN res(BN::BigInteger1());
    BN t = *this;

    do {
        if (power & 1)
            res = res * t;
        power >>= 1;
        if (power)
            t = t.qrt();
    } while (power);

    return res;
}

BN BN::powmod(uint64_t power, const BN& mod) const
{
    if (power == 0)
        return BN::BigInteger1();

    BN res(BN::BigInteger1());
    BN t = *this % mod;

    do {
        if (power & 1)
            res = res * t % mod;
        power >>= 1;
        if (power)
            t = t.qrt() % mod;
    } while (power);

    return res;
}

BN BN::powmod(const BN& power, const BN& mod) const {
    return expRightToLeft(power, mod);
}

BN BN::powmodBarrett(const BN& power, const BN& mod) const {
    if(power.isZero())
        return BN::BigInteger1();

    BN mu = reductionBarrettPrecomputation(mod);
    BN res(BN::BigInteger1());
    BN t = (*this) % mod;

    int len = power.bit_count();
    bt mask = 1;
    const bt *curr = &*power.data.begin();
    for(int i = 0; i < len; i++) {
        if(!mask) {
            mask = 1;
            ++curr;
        }
        if( (*curr) & mask)
            res = (res*t).reductionBarrett(mod, mu);

        if (i + 1 != len)
            t = t.qrt().reductionBarrett(mod, mu);
        mask <<= 1;
    }
    return res;
}

// Root operations

BN BN::sqrt() const
{
    if(isZero())
        return BN::BigInteger1();

    size_t rbc2 = (data.size() + 1) / 2 + 1;
    BN x(rbc2, 0);

    x.data.back() = 1;
    BN x0;

    do {
        x0 = x;
        x += *this / x;
        x >>= 1;
    } while(x0 > x);

    return x0;
}

BN BN::qrt() const
{
    if (data.size() < max_size_fast_mul)
        return fast_qrt();

    BN res(2 * data.size() + 1, 0);
    for (size_t i = 0; i < data.size(); ++i) {
        bt4 cuv = res.data[2 * i] + static_cast<bt2>(data[i]) * data[i];
        res.data[2 * i] = cuv;
        for (size_t j = i + 1; j < data.size(); ++j) {
            cuv = static_cast<bt4>(res.data[i + j]) +
                ((static_cast<bt4>(data[i]) * data[j]) << 1) +
                (cuv >> bz8);
            res.data[i + j] = cuv;
        }
        cuv = res.data[i + data.size()] + (cuv >> bz8);
        res.data[i + data.size()] = cuv;
        res.data[i + data.size() + 1] += (cuv >> bz8);
    }

    norm(res.data);
    return res;
}

BN BN::fast_qrt() const
{
    size_t n = data.size();

    BN result(n + n, 0);

    bt4 t = 0;
    for(size_t s = 0; s < n + n - 1; s++) {

        size_t start_index = s >= n ? s - n + 1 : 0;
        size_t end_index = min(n - 1, s);
        while (start_index < end_index) {
            bt2 m = static_cast<bt2>(data[start_index]) * data[end_index];
            t += m;
            t += m;
            ++start_index;
            --end_index;
        }
        if (start_index == end_index)
            t += static_cast<bt2>(data[start_index]) * data[end_index];

        result.data[s] = t;
        t = t >> bz8;
    }

    result.data[n + n - 1] = t;
    norm(result.data);
    return move(result);
}

// Left-Right rotations

BN BN::expRightToLeft(const BN& exponent, const BN& mod) const 
{
    if(exponent.isZero())
        return BN::BigInteger1();

    BN result(BN::BigInteger1());
    BN S = *this % mod;

    size_t len = exponent.bit_count();
    bt mask = 1;
    const bt *curr = &*exponent.data.begin();
    for(size_t i = 0; i < len; i++) {
        if(!mask) {
            mask = 1;
            ++curr;
        }
        if(*curr & mask)
            result = move(result * S % mod);

        if (i + 1 != len)
            S = move(S.qrt() % mod);
        mask <<= 1;
    }
    return result;
}

vector<BN> BN::expLeftToRightK_aryPrecomputation(const BN& mod) const 
{
    BN g = *this % mod;
    vector <BN> garr(KarySize);
    garr[0] = BN::BigInteger1();
    for(size_t i = 1; i < KarySize; i++) {
        garr[i] = garr[i-1] * g % mod;
    }
    return move(garr);
}

BN BN::expLeftToRightK_ary(const BN& exponent, const BN& mod, const vector<BN>& g) const 
{
    if(exponent.isZero())
        return BN::BigInteger1();

    BN A(BN::BigInteger1());
    for(int i = exponent.data.size() - 1; i >= 0; i--) {
        bt value = exponent.data[i];
        for (size_t b = bz - 1; b < bz; --b) {
            for(size_t k = 0; k < KaryBits; k++)
                A = A.qrt() % mod;
            A = A * g[(value >> KaryBits * b) & KaryMask] % mod;
        }
    }
    return A;
}

vector<BN> BN::expLeftToRightK_aryVarPrecomputation(const BN& mod, size_t K) const 
{
    int Kmax = (1 << K);
    BN g = *this % mod;

    vector <BN> garr(Kmax);
    garr[0] = BN::BigInteger1();

    for(int i = 1; i < Kmax; i++) {
        garr[i] = garr[i-1] * g % mod;
    }

    return move(garr);
}

BN BN::expLeftToRightK_aryVar(const BN& exponent, const BN& mod, const vector<BN>& g, size_t K) const 
{
    if(exponent.isZero())
        return BN::BigInteger1();

    BN A(BN::BigInteger1());

    int x = K;
    for(size_t i = exponent.data.size() * bz8 - 1; i >= K; i -= K) {
        x = i;
        for(size_t k = 0; k < K; k++)
            A = A.qrt() % mod;
        int curr = 0;
        for(size_t k = 0; k < K; k++) {
            curr <<= 1;
            curr |= exponent.get_bit(i-k);
        }
        A = A * g[curr] % mod;
    }

    uint32_t curr = 0;
    for(int i = x - K; i >= 0; i--) {
        A = A.qrt() % mod;
        curr <<= 1;
        curr |= exponent.get_bit(i);
    }
    return A * g[curr] % mod;
}

vector <BN> BN::expLeftToRightK_aryModifPrecomputation(const BN& mod) const 
{
    BN g = *this % mod;
    vector<BN> garr(KarySize);

    garr[0] = BN::BigInteger1();
    garr[1] = g; garr[2] = g.qrt() % mod;

    for(size_t i = 1; i < KarySize / 2; i++)
        garr[2 * i + 1] = garr[2 * i - 1] * garr[2] % mod;

    return move(garr);
}

BN BN::expLeftToRightK_aryMod(const BN& exponent, const BN& mod, const vector<BN>& g) const 
{
    if(exponent.isZero())
        return BN::BigInteger1();

    BN A(BN::BigInteger1());
    for(int i = exponent.data.size() - 1; i >= 0; i--) {
        for (size_t b = bz - 1; b < bz; --b) {
            bt ei = (exponent.data[i] >> KaryBits * b) & KaryMask;

            size_t hi = 0;
            if(ei != 0) {
                while(! (ei & 1)) {
                    ei >>= 1;
                    hi++;
                }
            }

            for(size_t k = 0; k + hi < KaryBits; k++)
                A = A.qrt() % mod;
            A = A * g[ei] % mod;
            for(size_t k = 0; k < hi; k++)
                A = A.qrt() % mod;
        }
    }

    return A;
}

vector<BN> BN::expSlidingWindowPrecomputation(const BN& mod, size_t k) const 
{
    size_t k_pow = 2 << (k-1);
    vector <BN> garr (k_pow);
    BN g = *this % mod;

    garr[0] = BN::BigInteger1(); 
    garr[1] = g;
    garr[2] = g.qrt() % mod;

    for(size_t i = 1; i < k_pow / 2; i++)
        garr[2 * i + 1] = garr[2 * i - 1] * garr[2] % mod;

    return move(garr);
}

// Sliding Window

BN BN::expSlidingWindow(const BN& exponent, const BN& mod, const vector<BN>& g, size_t K) const 
{
    BN A(BN::BigInteger1());
    int i = exponent.bit_count() - 1;

    while (i >= 0) {
        if(exponent.get_bit(i) == 0) {
            A = A.qrt() % mod;
            i--;
            continue;
        }
        int l = max(i - static_cast<int>(K) + 1, 0);
        while(exponent.get_bit(l) == 0)
            l++;

        int gx = 0;
        for(int j = i; j >= l; j--)
            gx = (gx << 1) | exponent.get_bit(j);
        for(int j = 0; j < i - l + 1; j++)
            A = A.qrt() % mod;
        A = A * g[gx] % mod;
        i = l - 1;
    }

    return A;
}

vector<BN> BN::expBest_SlidePrecomp(const BN& mod) const 
{
    vector <BN> garr (KarySize);
    BN mu = reductionBarrettPrecomputation(mod);
    BN g = this -> reductionBarrett(mod, mu);

    garr[0] = BN::BigInteger1();
    garr[1] = g;
    garr[2] = g.qrt().reductionBarrett(mod,mu);

    for(bt2 i = 1; i < KarySize/ 2; i++)
        garr[2 * i + 1] = (garr[2 * i - 1] * garr[2]).reductionBarrett(mod,mu);

    return move(garr);
}

BN BN::expBest_Slide(const BN& exponent, const BN& mod, const vector<BN>& g) const 
{
    BN A(BN::BigInteger1());
    BN mu = reductionBarrettPrecomputation(mod);
    int i = exponent.bit_count() - 1;
    int k = KaryBits;

    while (i >= 0) 
    {
        if(exponent.get_bit(i) == 0) {
            A = A.qrt().reductionBarrett(mod, mu);
            i--;
            continue;
        }

        int l = max(i - k + 1, 0);
        while(exponent.get_bit(l) == 0)
            l++;

        int gx = 0;
        for(int j = i; j >= l; j--)
            gx = (gx << 1) | exponent.get_bit(j);
        for(int j = 0; j < i - l + 1; j++)
            A = A.qrt().reductionBarrett(mod, mu);

        A = (A * (g[gx])).reductionBarrett(mod, mu);
        i = l - 1;
    }

    return A;
}

// Additional bit operations

size_t BN::count_zero_right() const noexcept
{
    if (data[0] & 1)
        return 0;
    if (isZero())
        return 0;

    size_t count = 0;
    while(!data[count])
        ++count;

    bt last = data[count];
    size_t result = count * bz8;
    while(!(last & 1)) {
        ++result;
        last >>= 1;
    }
    return result;
}

bool BN::get_bit(size_t index) const noexcept 
{
    if(index >= bz8 * data.size())
        return false;

    bt mask = 1;
    mask <<= (index % bz8);

    if (data[index / bz8] & mask)
        return true;
    return false;
}

uint64_t BN::get64() const noexcept 
{
    uint64_t result = 0;
    for(size_t i = min(data.size(), sizeof(uint64_t) / bz) - 1; i < sizeof(uint64_t); i--)
        result = (result << bz8) | data[i];
    return result;
}

const vector<bt> BN::raw() const noexcept 
{
    return data;
}

// Value check operations

bool BN::isZero() const noexcept
{
    if (data.size() > 1 || data[0])
        return false;
    return true;
}

bool BN::isEven() const noexcept
{
    if(data[0] & 1)
        return false;
    return true;
}

// Static default BigInteger instances 0 and 1.

const BN BN::BigInteger0() noexcept {
    static BN bn(0);
    return bn;
}

const BN BN::BigInteger1() noexcept {
    static BN bn(1);
    return bn;
}

// Random generator.
BN BN::random(size_t byteCount) 
{
    vector<bt> result(byteCount / bz);
    if (result.empty())
        result.emplace_back(rand() & bmax);
    else
        for (auto& i : result)
            i = rand() & bmax;

    return move(BN(move(result)));
}

// Print operations
string to_hexstring(const BN& bn) 
{
    string result;

    const auto& raw = bn.raw();
    for (auto i = raw.rbegin(); i != raw.rend(); ++i) {
        stringstream stream;
        stream << hex << setfill('0') << setw(bz * 2) << static_cast<uint32_t>(*i);
        string group = stream.str();
        result = result + group;
    }

    return result;
}

string to_string(BN bn) 
{
    stack<char> chars;

    do {
        chars.push(bn.modbase(10).get64() + '0');
        bn.divbaseappr(10);
    } while (!bn.isZero());

    string s;
    s.reserve(chars.size() + 1);
    while (!chars.empty()) {
        s = s + chars.top();
        chars.pop();
    }

    return s;
}

