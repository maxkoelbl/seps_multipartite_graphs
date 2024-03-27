"""
Proof of Proposition 3.2
"""

# The following functions provide the core algorithm

n = var('n')

def times2tp1(f):
    # Returns the generating function equivalent of multiplying its underlying polynomial by (2*k+1)
    return 2*t*f.derivative(t) + f

def recursion_checker(top, mid, bot):
    # Checks whether a recursive equation of the form
    # H1(k) = (2k+1) * a * H2(k) + alpha_0 * H3(k) + ... + alpha_1 * H3n(k)
    # holds using countable families of generating functions as inputs.
    # The families of generating functions should be of the form (n,t) --> HS_{n}(t).

    # Input:	top: generating function of a polynomial H1 of denominator degree d
    #			mid: generating function of a polynomial H2 of denominator degree d-1
    #			bot: list of the HS of the Ehrhart polynomials H30...H3n of denominator degree d-2
    # Output:	Either a list [a, alpha_0, ..., alpha_n] or None

    # This is based on an algorithm in Higashitani-Kummer-Michałek (2017).

    # Defining the variables t and
    t,a = var('t,alpha')
    alpha = [var('alpha'+str(i)) for i in range(len(bot))]
    test = a*times2tp1(mid)
    for j in range(len(bot)):
        test += alpha[j]*bot[j]
    test = (test/top).full_simplify()
    num = test.numerator()
    den = test.denominator()
    numCoeff = num.coefficients(t)
    denCoeff = den.coefficients(t)
    solveFor = [numCoeff[i][0]==denCoeff[i][0] for i in range(len(numCoeff))]

    solution = solve(solveFor,[a] + alpha)

    if len(solution) == 0:
        print("There are no solutions.")
        return None
    else:
        return solution[0]

# The following functions return the Ehrhart series of

def h(m, ring=QQ):
    # K_{m,n}
    gamma = lambda i: binomial(2*i,i)*binomial(m-1,i)*binomial(n-1,i)
    n,t = var('n'), var('t')
    deg = m+n-1
    hStar = 0
    for i in range(m):
        hStar += gamma(i)*t^i*(1+t)^(deg-2*i)
    return hStar/(1-t)^(deg+1)

def h_1(m, ring=QQ):
    # K_{1,m,n}
    gamma = lambda i: binomial(2*i,i)*binomial(m,i)*binomial(n,i)
    n,t = var('n'), var('t')
    deg = m+n
    hStar = 0
    for i in range(m+1):
        hStar += gamma(i)*t^i*(1+t)^(deg-2*i)
    return hStar/(1-t)^(deg+1)

def h_22(ring=QQ):
    # K_{2,2,n}
    gamma = lambda i: [1, 2*binomial(3*n+1,1), 2*binomial(3*n,2), 20*binomial(n,3)][i]
    n,t = var('n'), var('t')
    deg = 3+n
    hStar = 0
    for i in range(4):
        hStar += gamma(i)*t^i*(1+t)^(deg-2*i)
    return hStar/(1-t)^(deg+1)

def h_111(ring=QQ):
    #K_{1,1,1,n}
    gamma = lambda i: [1, 2*binomial(2*n+1,1), 6*binomial(n,2)][i]
    n,t = var('n,t')
    deg = 2+n
    hStar = 0
    for i in range(3):
        hStar += gamma(i)*t^i*(1+t)^(deg-2*i)
    return hStar/(1-t)^(deg+1)

print("(a) h{1,1,n} <- h{1,n}; [ h{1,n-1} ]")
print(
    recursion_checker(
        h_1(1)(n=n),
        h(1)(n=n),
        [
            h(1)(n=n-1)
        ]
    )
)
print("\n(b) h{1,1,n+1} <- h{1,1,n}; [ h{1,1,n-1} , h{1,n} ]")
print(
    recursion_checker(
        h_1(1)(n=n+1),
        h_1(1)(n=n),
        [
            h_1(1)(n=n-1),
            h(1)(n=n)
        ]
    )
)
print("\n(c) h{1,2,n} <- h{1,1,n}; [ h{1,1,n-1} , h{1,n} ]")
print(
    recursion_checker(
        h_1(2)(n=n),
        h_1(1)(n=n),
        [
            h_1(1)(n=n-1),
            h(1)(n=n)
        ]
    )
)
print("\n(d) h{1,2,n+1} <- h{1,2,n}; [ h{1,2,n-1} , h{1,1,n} , h{1,n+1} ]")
print(
    recursion_checker(
        h_1(2)(n=n+1),
        h_1(2)(n=n),
        [
            h_1(2)(n=n-1),
            h_1(1)(n=n),
            h(1)(n=n+1)
        ]
    )
)
print("\n(e) h{1,1,1,n} <- h{1,1,n}; [ h{1,1,n-1} , h{1,n}]")
print(
    recursion_checker(
        h_111()(n=n),
        h_1(1)(n=n),
        [
            h_1(1)(n=n-1),h(1)(n=n)
        ]
    )
)
print("\n(f) h{4,n} <- h{3,n}; [ h{3,n-1} , h{2,n} , h{1,n+1} ]")
print(
    recursion_checker(
        h(4)(n=n),
        h(3)(n=n),
        [
            h(3)(n=n-1),
            h(2)(n=n),
            h(1)(n=n+1)
        ]
    )
)
print("\n(g) h{3,n+1} <- h{3,n}; [ h{2,n} , (2k+1)h{2,n-1} , h{1,n+1} ]")
print(
    recursion_checker(
        h(3)(n=n+1),
        h(3)(n=n),
        [
            h(2)(n=n),
            times2tp1(h(2)(n=n-1)),
            h(1)(n=n+1)
        ]
    )
)
print("\n(h) h{2,2,n} <- h{1,2,n}; [ h{1,2,n-1} , h{1,1,n} , h{1,n+1} ]")
print(
    recursion_checker(
        h_22()(n=n),
        h_1(2)(n=n),
        [
            h_1(2)(n=n-1),
            h_1(1)(n=n),
            h(1)(n=n+1)
        ]
    )
)
print("\n(i) h_{1,3,n} <- h{1,2,n}; [h{1,2,n-1} , h{1,1,n} , h{1,n+1} ]")
print(
    recursion_checker(
        h_1(3)(n=n),
        h_1(2)(n=n),
        [
            h_1(2)(n=n-1),
            h_1(1)(n=n),
            h(1)(n=n+1)
        ]
    )
)
