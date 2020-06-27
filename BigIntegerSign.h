#ifndef _BIGINTEGERSIGN_H
#define _BIGINTEGERSIGH_H

#include "BigInteger.h"

#define BN BigInteger
#define BNsign BigIntegerSign

class BigIntegerSign
{
public:
    BN value;
    // false: positive
    // true: negative
    bool sign;

public:
    BNsign();
    BNsign(const BNsign& bn);
    BNsign(BNsign&& bn) noexcept;
    BNsign(const BN& bn, bool negative = false);
    BNsign(BN&& bn, bool negative = false) noexcept;
    BNsign& operator = (const BNsign&);
    BNsign& operator = (BNsign&&) noexcept;

    const BNsign operator + (const BNsign&) const;
    BNsign& operator += (const BNsign& );

    const BNsign operator - (const BNsign&) const;
    BNsign& operator -= (const BNsign& );

    const BNsign operator * (const BNsign&) const;
};

std::string to_string(const BNsign& bn);
std::string to_hexstring(const BNsign& bn);

#endif