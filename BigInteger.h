#ifndef _BIGINTEGER_H
#define _BIGINTEGER_H

#include <stdexcept>
#include <string>
#include <vector>

#define DOUBLE_BASE

#ifndef DOUBLE_BASE

using bt = uint8_t;
using bt2 = uint16_t;
using bt2s = int16_t;
using bt4 = uint32_t;
using bt4s = int32_t;
constexpr bt2 bsize = 256;
constexpr bt bmax = 255;

#else

using bt = uint16_t;
using bt2 = uint32_t;
using bt2s = int32_t;
using bt4 = uint64_t;
using bt4s = int64_t;
constexpr bt2 bsize = 65536;
constexpr bt bmax = 65535;

#endif

using namespace std;

constexpr bt bz = sizeof(bt);
constexpr bt bz8 = sizeof(bt) * 8;

enum class RadixBI {hex, dec};

struct divide_by_zero : logic_error {
    divide_by_zero();
};

class BigInteger {

    public:
        static const BigInteger BigInteger0() noexcept;
        static const BigInteger BigInteger1() noexcept;
        static BigInteger random(size_t byte);
    public:
        /*  CONSTRUCTORS */

        BigInteger();
        BigInteger(uint64_t base, int type);

        explicit BigInteger(uint64_t x);

        BigInteger(const BigInteger&);
        BigInteger(BigInteger&& big_int) noexcept;

        BigInteger(vector<bt>&&) noexcept;
        BigInteger(const vector<bt>&);

        BigInteger(vector<bt>&&, size_t rbc);
        BigInteger(const vector<bt>&, size_t rbc);

        BigInteger(const string &, RadixBI = RadixBI::hex);

        /*  CONSTRUCTORS END */

        // Swap two BigInteger instances.
        void swap(BigInteger& big_int) noexcept;

        // Elementary base operations
        BigInteger mulbt(size_t mul) const;
        BigInteger divbt(size_t div) const;
        BigInteger modbt(size_t mod) const;

        // Base multiplication
        const BigInteger mulbase(const bt&) const;
        BigInteger& mulbaseappr(const bt&);

        // Base division
        const BigInteger divbase(const bt&) const;
        BigInteger& divbaseappr(const bt&);

        // Base modulo
        const BigInteger modbase(const bt&) const;
        BigInteger& modbaseappr(const bt&);

        // Assignment operations
        BigInteger& operator = (const BigInteger&);
        BigInteger& operator = (BigInteger&&) noexcept;

        // Addition operations
        const BigInteger operator + (const BigInteger&) const;
        BigInteger& operator += (const BigInteger&);
        BigInteger& operator ++ ();

        // Subtraction operations
        const BigInteger operator - (const BigInteger&) const;
        BigInteger& operator -= (const BigInteger&);
        BigInteger& operator -- () noexcept;
        
        // Multiplication operations
        const BigInteger operator * (const BigInteger&) const;
        // naive multiplication O(n ^ 2)
        const BigInteger naive_mul(const BigInteger&) const;
        // quick multiplication O(n ^ 2)
        const BigInteger fast_mul(const BigInteger&) const;
        // Karatsuba Algorithm O(n ^ lg3))
        const BigInteger karatsuba_mul(const BigInteger&) const;

        // Division operations
        void divmod(const BigInteger&, BigInteger& div, BigInteger& mod) const;
        const BigInteger operator / (const BigInteger&) const;
        const BigInteger operator % (const BigInteger&) const;

        // Right shift operations
        const BigInteger operator >> (size_t shift) const;
        BigInteger& operator >>= (size_t shift);

        // Left shift operations
        const BigInteger operator << (size_t shift) const;
        BigInteger& operator <<= (size_t shift);

        // Conditional operations
        bool operator <  (const BigInteger&) const noexcept;
        bool operator <= (const BigInteger&) const noexcept;
        bool operator >  (const BigInteger&) const noexcept;
        bool operator >= (const BigInteger&) const noexcept;
        bool operator == (const BigInteger&) const noexcept;
        bool operator != (const BigInteger&) const noexcept;
        bt   operator [] (size_t index_base) const noexcept;

        // Digit information
        size_t digit_count() const noexcept;
        size_t bit_count()   const noexcept;

        // Power Operations (Barret's reduction algorithm)
        BigInteger reductionBarrett(const BigInteger& mod,const BigInteger& mu) const;
        BigInteger pow (uint64_t) const;
        BigInteger powmod (uint64_t power, const BigInteger& mod) const;
        BigInteger powmod (const BigInteger& power, const BigInteger& mod) const;
        BigInteger powmodBarrett(const BigInteger& power, const BigInteger& mod) const;

        // Root operations
        BigInteger sqrt()       const;
        BigInteger qrt()        const;
        BigInteger fast_qrt()   const;

        // Value check operations
        bool isZero() const noexcept;
        bool isEven() const noexcept;

        BigInteger expRightToLeft(const BigInteger& power, const BigInteger& mod) const;

        vector<BigInteger> expLeftToRightK_aryPrecomputation(const BigInteger& mod) const;
        BigInteger expLeftToRightK_ary(const BigInteger& exponent, const BigInteger& mod, const vector<BigInteger>& g) const;

        vector<BigInteger> expLeftToRightK_aryVarPrecomputation(const BigInteger& mod, size_t k) const;
        BigInteger expLeftToRightK_aryVar(const BigInteger&, const BigInteger&, const vector<BigInteger>&, size_t k) const;

        vector<BigInteger> expLeftToRightK_aryModifPrecomputation(const BigInteger&) const;
        BigInteger expLeftToRightK_aryMod(const BigInteger&, const BigInteger&, const vector<BigInteger>&) const;

        vector<BigInteger> expSlidingWindowPrecomputation(const BigInteger&, size_t k) const;
        BigInteger expSlidingWindow(const BigInteger&, const BigInteger&, const vector<BigInteger>&, size_t k) const;

        //for best result:
        vector<BigInteger> expBest_SlidePrecomp(const BigInteger& mod) const;
        BigInteger expBest_Slide(const BigInteger& exponent, const BigInteger& mod, const vector<BigInteger>& g) const;

        const vector<bt> raw() const noexcept;

    private:
        vector<bt> data;
};

string to_string(BigInteger);
string to_hexstring(const BigInteger &);

#endif