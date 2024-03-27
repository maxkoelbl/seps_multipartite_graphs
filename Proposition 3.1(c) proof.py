"""
Proof of Proposition 3.1(c)

We use the formula from Theorem 2.5 to compute h*_{2,2,n}.
Then we simplify it and compare it with the reference formula from Proposition 3.1(c) coefficient-wise.
"""

# Defining variable names.

n,t,nu = var('n,t,nu')
a,b,c = var('a,b,c')
x,y = var('x,y')
j2,j3,j4,j5,j6 = var('j2,j3,j4,j5,j6')

# The function p as given in Theorem 2.5.

p(x,y,i,j) = binomial(x-y-1,i)*binomial(y-1,j)*binomial(y+i-j-1,i)

# The function c with four explicit numbers of inputs following the definition in Section 2.2.

c_3(c1,c2,c3) = sum(
    binomial(c1+j2-1,j2)
    *binomial(c2-j2+c3-2,c3-1),
    j2, 0, c2-1
)

c_4(c1,c2,c3,c4) = sum(
    sum(
        binomial(c1+j2-1,j2)
        *binomial(c3-j3+c4-2,c4-1)
        *binomial(c2-j2+j3-1,j3)
        , j3, 0, c3-1
    ), j2, 0, c2-1
)

c_5(c1,c2,c3,c4,c5) = sum(
    sum(
        sum(
            binomial(c1+j2-1,j2)
            *binomial(c4-j4+c5-2,c5-1)
            *binomial(c2-j2+j3-1,j3)
            *binomial(c3-j3+j4-1,j4),
            j4, 0, c4-1
        ), j3, 0, c3-1
    ), j2, 0, c2-1
)

c_6(c1,c2,c3,c4,c5,c6) = sum(
    sum(
        sum(
            sum(
                binomial(c1+j2-1,j2)
                *binomial(c5-j5+c6-2,c6-1)
                *binomial(c2-j2+j3-1,j3)
                *binomial(c3-j3+j4-1,j4)
                *binomial(c4-j4+j5-1,j5),
                j5, 0, c5-1
            ), j4, 0, c4-1
        ), j3, 0, c3-1
    ), j2, 0, c2-1
)

# The function h^(i)_{a,b,c} as given in Theorem 2.5.

h_i(a,b,c,t) = sum(
    sum(
        p(a+b+c,a,i,j)
        *binomial(b+c+j-i-2,j-1)
        *(t^(i+j+1) + t^(a+b+c-i-j-2)),
        j,1,a-1
    ), i, 0, b+c-1
) + sum(
    sum(
        p(a+b+c,b,i,j)
        *binomial(a+c-i+j-2,a+c-i-1)
        *(t^(i+j) + t^(a+b+c-i-j-1)),
        j,1,b-1
    ), i, 0, a+c-1
) + sum(
    sum(
        p(a+b+c,c,i,j)
        *binomial(a+b-i+j-2,a+b-i-1)
        *(t^(i+j) + t^(a+b+c-i-j-1)),
        j,1,c-1
    ), i, 0, a+b-1
)

# Substituting a=b=2 and c=n.

h_i_22n(n,t) = h_i(2,2,n,t)

# All 13 cases of q(nu_1,nu_2,nu_3) r(nu_1,nu_2,nu_3) (t^{nu_1+nu_2+nu_3} + t^{a+b+c-nu_1-nu_2-nu_3-1}).
# Since a=b=2, we can insert 1 for nu_1, a-nu_1, nu_2, and b-nu_2 whenever they appear.

# nu_1 = 0 , nu_2 = 2 , nu_3 = n

part1(n,t) = c_3(2,2,n) * (t^(2+n) + t)

# nu_1 = 0 , nu_2 = 2 , nu_3 = 0

part2(n,t) = c_3(2,2,n) * (t^2 + t^(1+n))

# nu_1 = 0 , nu_2 = 0 , nu_3 = n

part3(n,t) = c_3(2,n,2) * (t^n + t^3)

# nu_1 = 0 , nu_2 = 2 , nu_3 not in {0,n}

part4(nu,n,t) = binomial(n,nu) * c_4(nu,2,2,n-nu) * (t^(2+nu) + t^(1-nu+n))

# nu_1 = 0 , nu_2 = 1 , nu_3 = n

part5(n,t) = 2 * c_4(1,2,n,1) * (t^(1+n) + t^2)

# nu_1 = 1 , nu_2 = 2 , nu_3 = 0

part6(n,t) = c_4(1,n,2,1) * (t^3 + t^n)

# nu_1 = 1 , nu_2 = 0 , nu_3 = n

part7(n,t) = c_4(1,2,n,1) * (t^(1+n) + t^2)

# nu_1 = 0 , nu_2 = 1 , nu_3 not in {0,n}

part8(nu,n,t) = 2 * binomial(n,nu) * c_5(1,nu,2,1,n-nu) * (t^(1+nu) + t^(2-nu+n))

# nu_1 = 1 , nu_2 = 0 , nu_3 not in {0,n}

part9(nu,n,t) = binomial(n,nu) * c_5(1,nu,2,1,n-nu) * (t^(1+nu) + t^(2-nu+n))

# nu_1 = 1 , nu_2 = 2 , nu_3 not in {0,n}

part10(nu,n,t) = binomial(n,nu) * c_5(1,n-nu,2,1,nu) * (t^(3+nu) + t^(-nu+n))

# nu_1 = 1 , nu_2 = 1 , nu_3 = 0

part11(n,t) = 2 * c_5(1,1,n,1,1) * (t^2 + t^(1+n))

# nu_1 = 1 , nu_2 = 1 , nu_3 = n

part12(n,t) = 2 * c_5(1,1,n,1,1) * (t^(2+n) + t)

# nu_1 = 1 , nu_2 = 1 , nu_3 not in {0,n}

part13(nu,n,t) = 2 * binomial(n,nu) * (
    c_6(1,1,n-nu,1,1,nu) +
    c_6(1,n-nu,1,1,nu,1) +
    c_6(1,1,nu,1,1,n-nu)
) * (t^(2+nu) + t^(1-nu+n))

# Summing up.

h_ii_22n(n,t) = (
    part1(n,t) +
    part2(n,t) +
    part3(n,t) +
    part5(n,t) +
    part6(n,t) +
    part7(n,t) +
    part11(n,t) +
    part12(n,t)
) + sum(
    part4(nu,n,t) +
    part8(nu,n,t) +
    part9(nu,n,t) +
    part10(nu,n,t) +
    part13(nu,n,t),
    nu, 1, n-1
)

# Assembling the h*-polynomial.

hstar(n,t) = (h_i_22n(n,t) + h_ii_22n(n,t)).full_simplify()

# Two new variables

i,j = var('i,j')

# The following functions are pieces we will use to assemble hstar in an alternative way.
# We will run an automatic check to verify the second definition is indeed identical

f_1(n) = -1/6*(2*(n - 1)^3 - 3*(n - 1)^2*n + 9*(n - 1)^2 - 3*(n - 1)*n - 12*n^2 + n - 7)
f_2(n) = 1/2*(3*(n - 1)^2 + 21*n + 1)
f_3(n) = 3*n + 2
f_4(n,i) = binomial(n + 1, i)
f_5(n,j) = 3*j*binomial(-j + n + 1, 2)*binomial(n - 1, j)
f_6(n,j) = binomial(j + 2, 3)*binomial(n - 1, j)
f_7(n,j) = binomial(-j + n + 2, 3)*binomial(n - 1, j)
f_8(n,j) = -3*(j*binomial(j + 1, 2)*binomial(n - 1, j) - n*binomial(j + 1, 2)*binomial(n - 1, j))
f_9(n,nu) = -3/2*((nu - 1)^2 - 2*nu^2 - nu - 1)*binomial(n, nu)
f_10(n,nu) = (3*n*nu*binomial(n, nu) - 3*nu^2*binomial(n, nu) + (3*n + 2)*binomial(n, nu))
f_11(n,nu) = -1/2*((n - nu - 1)^2 - 2*(n - nu)*n + 2*(n - nu)*nu - n + nu - 1)*binomial(n, nu)

# Assembling hstar.

hstar_alt(n,t) = (
    f_1(n)*(t^n + t^3)
    + f_2(n)*(t^(n+1) + t^2)
    + f_3(n)*(t^(n+2) + t)
    + sum(
        f_4(n,i) * (t^(i + 2) + t^(n + 1 - i)),
        i, 0, n + 1
    )
    + sum(
        f_4(n,i) * (t^(i + 1) + t^(n + 2 - i)),
        i, 0, n + 1
    )
    + sum(
        f_5(n,j) * (t^(j + 2) + t^(n + 1 - j)),
        j, 1, n - 1
    )
    + sum(
        f_6(n,j) * (t^(j) + t^(n + 3 - j)),
        j, 1, n - 1
    )
    + sum(
        f_7(n,j) * (t^(j + 3) + t^(n - j)),
        j, 1, n - 1
    )
    + sum(
        f_8(n,j) * (t^(j + 1) + t^(n + 2 - j)),
        j, 1, n - 1
    )
    + sum(
        f_9(n,nu) * (t^(nu + 1) + t^(n - nu + 2))
        + f_10(n,nu) * (t^(nu + 2) + t^(n - nu + 1))
        + f_11(n,nu) * (t^(nu + 3) + t^(n - nu)),
        nu, 1, n - 1
    )
)

# Checking the for correctness.

print( "hstar_alt is indeed equivalent to hstar:" , bool(hstar(n,t) == hstar_alt(n,t)) )

# Based on hstar_alt, we define the initialise the coefficients of hstar.
# We start with a number of useful intermediary functions.

g_1(n) = f_4(n,n+1)
g_2(n) = f_3(n) + f_4(n,0) + f_4(n,n) + f_4(n,n+1)
g_3(n) = f_2(n) + f_4(n,0) + f_4(n,n)
g_4(n) = f_1(n)
g_5(n,i) = f_6(n,i)
g_6(n,i) = f_4(n,i) + f_8(n,i) + f_9(n,i)
g_7(n,i) = f_4(n,i) + f_5(n,i) + f_10(n,i)
g_8(n,i) = f_7(n,i) + f_11(n,i)

# The simplified h*-polynomial:

hstar_simplified(n,t) = (
      g_1(n) * (t^0 + t^(n + 3))
    + g_2(n) * (t^1 + t^(n + 2))
    + g_3(n) * (t^2 + t^(n + 1))
    + g_4(n) * (t^3 + t^(n + 0))

    + sum(
          g_5(n,i) * (t^(i + 0) + t^(n - i + 3))
        + g_6(n,i) * (t^(i + 1) + t^(n - i + 2))
        + g_7(n,i) * (t^(i + 2) + t^(n - i + 1))
        + g_8(n,i) * (t^(i + 3) + t^(n - i + 0)),
        i, 1, n - 1
    )
)

# Now we can initialise the coefficients hstar_0, hstar_1, hstar_2, h_star_3, and hstar_i for 3 < i ≤ ceil((n+3)/2).
# This assumes that n ≥ 4. We will check 1 ≤ n ≤ 3 explicitly afterwards.

hstar_0(n) = g_1(n)
hstar_1(n) = g_2(n) + g_5(n,1) + g_8(n,n-1)
hstar_2(n) = g_3(n) + g_5(n,2) + g_6(n,1) + g_8(n,n-2) + g_7(n,n-1)
hstar_3(n) = g_4(n) + g_5(n,3) + g_6(n,2) + g_7(n,1) + g_8(n,n-3) + g_7(n,n-2) + g_6(n,n-1)
hstar_i(n,i) = g_5(n,i) + g_6(n,i-1) + g_7(n,i-2) + g_8(n,i-3) + g_8(n,n-i) + g_7(n,n-i+1) + g_6(n,n-i+2) + g_5(n,n-i+3)

# Now we check if these coefficients match with the reference formula given in the proposition.

factor_1(n) = 2*binomial(3*n+1,1)
factor_2(n) = 2*binomial(3*n,2)
factor_3(n) = 20*binomial(n,3)

ref_0(n) = binomial(n + 3,0)
ref_1(n) = binomial(n + 3,1) + factor_1(n)*binomial(n + 1, 0)
ref_2(n) = binomial(n + 3,2) + factor_1(n)*binomial(n + 1, 1) + factor_2(n)*binomial(n - 1, 0)
ref_3(n) = binomial(n + 3,3) + factor_1(n)*binomial(n + 1, 2) + factor_2(n)*binomial(n - 1, 1) + factor_3(n)*binomial(n - 3, 0)
ref_i(n,i) = binomial(n + 3,i) + factor_1(n)*binomial(n + 1, i-1) + factor_2(n)*binomial(n - 1, i-2) + factor_3(n)*binomial(n - 3, i-3)

print( "0-th coefficients match:" , bool(ref_0(n) == hstar_0(n)) )
print( "1-st coefficients match:" , bool(ref_1(n) == hstar_1(n)) )
print( "2-nd coefficients match:" , bool(ref_2(n) == hstar_2(n)) )
print( "3-rd coefficients match:" , bool(ref_3(n) == hstar_3(n)) )
print( "i-th coefficients match:" , bool(ref_i(n,i) == hstar_i(n,i)) )

# Lastly, we check the cases 1 ≤ n ≤ 3 explicitly.

reference(n,t) = factor_3(n)*(1+t)^(n-3)*t^3 + factor_2(n)*(1+t)^(n-1)*t^2 + factor_1(n)*(1+t)^(n+1)*t + (1+t)^(n+3)

print( "Equivalence holds for n=1:" , bool(reference(1,t) == hstar_simplified(1,t)) )
print( "Equivalence holds for n=2:" , bool(reference(2,t) == hstar_simplified(2,t)) )
print( "Equivalence holds for n=3:" , bool(reference(3,t) == hstar_simplified(3,t)) )
