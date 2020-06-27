#ifndef _GCD_
#define _GCD_

#include <algorithm>
#include <iostream>
#include <vector>

#include "gcd.h"

BN gcdEuclidean(BN a,BN b)
{
    while(!b.isZero()) {
        a = a % b;
        a.swap(b);
    }
    return a;
}

BN gcdBinary(BN a,BN b)
{
    if(a.isZero())
        return b;

    if(b.isZero())
        return a;

    size_t acount = a.count_zero_right();
    size_t bcount = b.count_zero_right();

    a >>= acount;
    b >>= bcount;

    if(a < b)
        swap(a, b);

    while(!a.isZero()) {
        if(a >= b) {
            a -= b;
            a >>= a.count_zero_right();
        } else {
            b -= a;
            b >>= b.count_zero_right();
        }
    }
    return move(b <<= min(acount,bcount));
}

BN gcdInverseEuclidean(const BN& a, const BN& mod)
{
    BN startMod(mod);
    BN A(a % mod), B(mod);

    if(A.isZero())
        return move(BN::BigInteger0());

    BNsign x0;
    BNsign x1(BN::BigInteger0());
    BNsign x2(BN::BigInteger1());
    BN Div, Mod;
    while (true) {
        B.divmod(A, Div, Mod);
        if (Mod.isZero())
            break;
        x0 = move(x1);
        x1 = move(x2);

        x2 = x0 - x1 * (BNsign)Div;

        B = move(A);
        A = move(Mod);
    }

    if(A != BN::BigInteger1())
        return move(BN::BigInteger0());
    if(x2.sign && !x2.value.isZero())
        return move(startMod -= x2.value);

    return move(x2.value);
}

BN gcdInverseEuclideanBinary(BN xx, BN mod)
{
    BN& yy = mod;
    BN g(BN::BigInteger1());
    size_t xcount=xx.count_zero_right();
    size_t ycount=yy.count_zero_right();
    if (xcount && ycount)
        return BN::BigInteger0();

    xx >>= min(xcount,ycount);
    yy >>= min(xcount,ycount);

    BN u=xx;
    BN v=yy;
    BNsign x=xx;
    BNsign y=yy;

    BNsign a(BN::BigInteger1());
    BNsign b(BN::BigInteger0());
    BNsign c(BN::BigInteger0());
    BNsign d(BN::BigInteger1());

    do {
        size_t uZeros = u.count_zero_right();
        size_t vZeros = v.count_zero_right();
        for (size_t i = 0; i < uZeros; ++i) {
            if(!a.value.isEven() || !b.value.isEven()) {
                a += y;
                b -= x;
            }
            a.value >>= 1;
            b.value >>= 1;
        }
        for (size_t i = 0; i < vZeros; ++i) {
            if(!c.value.isEven() || !d.value.isEven()) {
                c += y;
                d -= x;
            }
            c.value >>= 1;
            d.value >>= 1;
        }
        u >>= uZeros;
        v >>= vZeros;
        if(u >= v) {
            u -= v;
            a -= c;
            b -= d;
        } else {
            v -= u;
            c -= a;
            d -= b;
        }
    }
    while(!u.isZero());

    if(v != BN::BigInteger1())
        return move(BN::BigInteger1());

    if(c.sign && !c.value.isZero())
        return move(mod -= c.value % mod);

    return move(c.value % mod);
}

vector <BN> multi_inverse(const vector <BN> &x, const BN &mod)
{
    int count=x.size();
    vector <BN> a;
    vector <BN> x_inverse;

    a.reserve(count);
    x_inverse.reserve(count);

    a.push_back(x[0]);

    for( vector<BN>::const_iterator xi=x.begin()+1;
        xi!=x.end();
        xi++)
        a.push_back(a.back() * (*xi) % mod);

    BN curr = gcdInverseEuclideanBinary(a[0],mod);
    for(int i=x.size()-1;i>0;i--)
    {
        x_inverse.push_back(a[i-1] * curr % mod);   
        curr = x[i] * curr % mod;
    }
    x_inverse.push_back(curr);               
    reverse(x_inverse.begin(),x_inverse.end());

    return x_inverse;
}

BN Garner(const std::vector <BN>& m, const std::vector <BN>& v)
{
    if(m.size()!=v.size())
        throw invalid_argument("Garner: Sizes of arrays are different");

    int t = m.size();
    vector<BN> C(t);
    BN u;
    for(int i=1;i<t;i++)
    {
        C[i] = BN::BigInteger1();
        for(int j=0; j<i; j++)
        {
            u = gcdInverseEuclideanBinary(m[j],m[i]);
            C[i] = u * C[i] % m[i];
        }
    }
    u = v[0];
    BN x = u;
    BN Pm(BN::BigInteger1());
    for(int i=1;i<t;i++)
    {
        BN xmi = x%m[i];
        BN vix = (v[i] >= xmi ? v[i] - xmi : m[i] + v[i] - xmi);
        u = vix * C[i] % m[i];
        Pm = Pm * m[i-1];
        x += u * Pm;
    }
    return x;
}

BN CTO(const std::vector <BN>& m, const std::vector <BN>& v)
{
    if(m.size()!=v.size())
        throw invalid_argument("CTO: Sizes of arrays are different");

    int t = m.size();
    BN M = m[0];
    for (int i=1;i<t;i++)
        M = M * m[i];
    BN x(BN::BigInteger0());
    for (int i=0;i<t;i++)
    {
        BN Mi = M/m[i];
        x += v[i] * (M/m[i]) * gcdInverseEuclideanBinary(M/m[i],m[i]) % M;
    }
    return x % M;
}

bool isPrimeMillerRabin(const BN& n, size_t k) 
{
    if (n.isEven()) return false;

    BN n1(n);
    --n1;
    size_t s = n1.count_zero_right();
    BN t = n1 >> s;
    
    for (size_t i = 0; i < k; ++i) 
    {
        BN a(uint64_t(i + 2));
        BN x = a.powmod(t, n);
        if (x == BN::BigInteger1() || x == n1)
            continue;
        
        bool stop = false;
        for (size_t j = 0; !stop && j < s - 1; ++j) 
        {
            x = x.qrt() % n;
            if (x == BN::BigInteger1())
                return false;
            if (x == n1)
                stop = true;
        }

        if (!stop)
            return false;
    }
    return true;
}


#endif