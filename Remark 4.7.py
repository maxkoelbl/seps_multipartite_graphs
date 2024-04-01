"""
On Remark 4.7

We run the formula from Corollary 4.5 for 1 ≤ m ≤ 10
"""

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

# The function that returns the Ehrhart series of the symmetric edge polytopes of K_{a,b}

def h(m, ring=QQ):
    # K_{m,n}
    gamma = lambda i: binomial(2*i,i)*binomial(m-1,i)*binomial(n-1,i)
    n,t = var('n'), var('t')
    deg = m+n-1
    hStar = 0
    for i in range(m):
        hStar += gamma(i)*t^i*(1+t)^(deg-2*i)
    return hStar/(1-t)^(deg+1)

# The loop that shows the alpha_i.

for m in [1..10]:
    print("\nm = ",m)
    for parameter in recursion_checker(
            h(m+1)(n=n+1), 
            h(m)(n=n+1),
            [
                h(m-i)(n=n+i) for i in [0..m-1]
            ]
        ):
        print(parameter)
