
struct PolynomialModP
    poly::Polynomial
    prime::Integer
    PolynomialModP() = new(Polynomial())
    PolynomialModP(prime::Integer) = new(Polynomial(),prime)
    PolynomialModP(poly::Polynomial,prime::Integer,dummy::Int) = new(poly,prime)
end

"""
Construct polynomial with a single term or multiple terms.
"""
function PolynomialModP(poly::Polynomial,prime::Integer)
    terms = Term[]
    for t in poly.terms
        if mod(t.coeff,prime) != 0
            push!(terms, mod(t,prime))
        end
    end
    p = Polynomial(terms,1)
    return PolynomialModP(p,prime,1)
end
"""
Push a new term into the polynomial.
"""
#Note that ideally this would throw and error if pushing another term of degree that is already in the polynomial
function push!(p::PolynomialModP, t::Term) 
    iszero(t) && return #don't push a zero
    push!(p.poly.terms,t)
end

"""
Pop the leading term out of the polynomial.
"""
pop!(p::PolynomialModP)= popat!(p.poly.terms,1)
"""
Check if the polynomial is zero.
"""
iszero(p::PolynomialModP)::Bool = isempty(p.poly.terms)

"""
The negative of a polynomial.
"""

function -(p::PolynomialModP)::PolynomialModP
    pol = p.poly
    pol_out = Polynomial()
    for i in pol
        push!(pol_out,Term(-i.coeff,i.degree))
    end
    return PolynomialModP(pol_out,p.prime)
end
"""
Leading term of PolynomialModP
"""
function leading(p::PolynomialModP)::Term 
    return p.poly.terms[1]
end
"""
Degree of PolynomialModP
"""
function degree(p::PolynomialModP)::Int
    if iszero(p)
        return 0
    else
        return leading(p.poly).degree
    end
end




"""
Add two polynomials.
"""
function +(p1::PolynomialModP, p2::PolynomialModP)::PolynomialModP
    prime = p1.prime
    p3 = PolynomialModP(prime)
    while !iszero(p1) && !iszero(p2)
        t1, t2 = leading(p1), leading(p2) 
        if t1.degree == t2.degree
            degree = t1.degree
            lead_p1 = pop!(p1)
            lead_p2 = pop!(p2)
            if mod(t1.coeff+t2.coeff,prime) != 0
                push!(p3, Term(mod(lead_p1.coeff+lead_p2.coeff,prime),lead_p1.degree))
            end
        elseif t1.degree < t2.degree
            push!(p3,pop!(p2))
        else
            push!(p3,pop!(p1))
        end
    end
    while !iszero(p1)
        push!(p3,pop!(p1))
    end
    while !iszero(p2)
        push!(p3,pop!(p2))
    end
    return p3
end

"""
Add a polynomial and a term.
"""
+(p::PolynomialModP, t::Term) = PolynomialModP(p.polynomial + Polynomial(t),p.prime)
+(t::Term, p::PolynomialModP) = PolynomialModP(p.polynomial + t,p.prime)

"""
Add a polynomial and an integer.
"""
+(p::PolynomialModP, n::Int) = PolynomialModP(p.polynomial + Term(n,0),p.prime)
+(n::Int, p::PolynomialModP) = PolynomialModP(p.polynomial + Term(n,0),p.prime)
"""
Subtraction of two polynomials.
"""

function -(p1::PolynomialModP,p2::PolynomialModP)::PolynomialModP
    @assert p1.prime == p2.prime
    pol1 = p1.poly
    pol2 = p2.poly
    neg_pol2 = -pol2
    pol_out = pol1 + neg_pol2
    return PolynomialModP(pol_out,p1.prime)
end

"""
Multiply two polynomials.
"""
function *(p1::PolynomialModP, p2::PolynomialModP)::PolynomialModP
    @assert p1.prime == p2.prime
    p_out = Polynomial()
    for t in p1.poly
        p_out = p_out + (t * p2.poly)
    end
    return PolynomialModP(p_out,p1.prime)
end

"""
Division of 2 PolynomialModP types
"""
function divide(sp1::PolynomialModP,sp2::PolynomialModP)
    @assert sp1.prime == sp2.prime
    iszero(sp2.poly) && throw(DivideError())
    p1, p2 = deepcopy(sp1), deepcopy(sp2)
    q = PolynomialModP(sp1.prime)
    prev_degree = degree(p1)
    while degree(p1) ≥ degree(p2)
        h_sub = (leading(p1)÷leading(p2))(sp1.prime)
        h_sub_2 = Polynomial(h_sub)
        h = PolynomialModP(h_sub_2,sp1.prime )  #syzergy 
        p1 = p1 - h*p2
        q = q + h
        prev_degree == degree(p1) && break
        prev_degree = degree(p1)
    end
    p1
    return q, p1
end

"""
Show function for PolynomialModP
"""
function show(io::IO, p::PolynomialModP) 
    p = deepcopy(p.poly)
    if iszero(p)
        print(io,"0")
    else
        n = length(p.terms)
        for (i,t) in enumerate(p.terms)
            print(io, t, i != n ? " + " : "")
        end
    end
end
"""
Derivative function for PolynomialModP
"""
function derivative(p::PolynomialModP)::PolynomialModP
    der_p = Polynomial()
    for term in p.poly
        der_term = derivative(term)
        !iszero(der_term) && push!(der_p,der_term)
    end
    return PolynomialModP(der_p,p.prime)
end
"""
Creates the unit polynomial.
"""
#one(prime::Integer)::PolynomialModP = PolynomialModP(Polynomial(one(Term())),prime)
#one(p::PolynomialModP) = one(typeof(p))

"""
Power of a polynomial mod prime.
"""
function pow_mod(p::PolynomialModP, n::Int)
    n < 0 && error("No negative power")
    out = deepcopy(p)
    if n == 1
        return out
    else 
        for _ in 2:n
            out *= p
        end
        return out
    end
end

"""
The extended euclid algorithm for polynomials modulo prime.
"""
#Trying to divide by 1 (r=1), need to fix and find out how to go from r being 1 to becoming 0. 
#Maybe a remainder of 1 (0*prime + 1), so then you do euclid(1,0).

function extended_euclid_alg(a::PolynomialModP, b::PolynomialModP)
    @assert a.prime == b.prime
    if a.poly == Polynomial(Term(0,0))
        return b, 0, 1
    elseif b.poly == Polynomial(Term(0,0))
        return a, 1, 0
    else
        old_r, r = a, b
        old_s, s = PolynomialModP(Polynomial(Term(1,0)),a.prime), PolynomialModP(Polynomial(Term(0,0)),a.prime)
        old_t, t = PolynomialModP(Polynomial(Term(0,0)),a.prime), PolynomialModP(Polynomial(Term(1,0)),a.prime)
        while !iszero(r)
            q = divide(old_r, r) |> first
            old_r, r = r, old_r-q*r
            old_s, s = s, old_s-q*s
            old_t, t = t, old_t-q*t
        end
        g, s, t = old_r, old_s, old_t
        return g, s, t  
    end
end

function test()
    a = PolynomialModP(Polynomial([Term(4,3),Term(2,2),Term(3,1),Term(2,0)]),7)
    b = PolynomialModP(Polynomial([Term(6,8),Term(3,4),Term(5,2),Term(2,1)]),7)
    ans = PolynomialModP(Polynomial([Term(3,11),Term(5,10),Term(4,9),Term(5,8),Term(5,7),Term(6,6),Term(1,5),Term(3,4),Term(5,3),Term(2,2),Term(4,1)]),7)
    @time mult_ans = a*b
    if ans.poly == mult_ans.poly
        return "Multiplication works"
    else
        return "Multiplication doesnt work"
    end
end
"""
Chinese remainder theorem on integers
"""


function chineseremainder(n::Array, a::Array)
    Π = prod(n)
    return mod(sum(ai * invmod(Π ÷ ni, ni) * (Π ÷ ni) for (ni, ai) in zip(n, a)), Π)
end


function CRT(a::Vector{PolynomialModP},n::Vector{Int64})::Polynomial
    c = Polynomial(Term(0,0))
    for k in max(degree(a[1]),degree(a[2])):-1:0
        if k == degree(a[1])
            if iszero(a[1])
                ak = 0
            else
                ak = pop!(a[1]).coeff 
            end
        else
            ak = 0
        end
        if k == degree(a[2])
            if iszero(a[2])
                bk = 0
            else
                bk = pop!(a[2]).coeff 
            end
        else
            bk = 0
        end
        ck = chineseremainder(n,[ak,bk])
        push!(c.terms,Term(ck,k))
    end
    return c
end

function *(a::Polynomial,b::Polynomial)
    a_height = 0
    for i in a.terms
        if abs(i.coeff) > a_height
            a_height = abs(i.coeff)
        end
    end
    b_height = 0
    for j in b.terms
        if abs(j.coeff) > b_height
            b_height = abs(j.coeff)
        end
    end
    B = 2*a_height*b_height*min(degree(a)+1,degree(b)+1)
    M = 3
    ap = PolynomialModP(a,M)
    bp = PolynomialModP(b,M)
    c = ap*bp
    while M < B
        p = nextprime(M+1)
        ap = PolynomialModP(a,p)
        bp = PolynomialModP(b,p)
        c2 = ap*bp
        c = PolynomialModP(CRT([c,c2],[M,p]),p)
        M = M*p
    end
    return smod(c)
end

"""
Symmetric mod to allow for negative terms
"""
function smod(a::PolynomialModP)::PolynomialModP
    prime = a.prime 
    out = PolynomialModP(a.prime)
    for i in a.poly
        if mod(i.coeff,prime) > prime//2
            i = Term(i.coeff - prime, i.degree)
            push!(out.poly,i)
        else 
            push!(out,i)
        end
    end
    return out 
end


function pow_mod(a::PolynomialModP, n::Int)::PolynomialModP
    n_bit = bitstring(n) 
    n_bit_split = split(n_bit,"")
    current = a
    out = PolynomialModP(Polynomial(Term(1,0)),a.prime)
    highest_one = 0
    for i in 1:64
        if n_bit_split[i] == "1"
            highest_one = i
            break
        end
    end
    for j in 64:-1:highest_one
        if n_bit_split[j] == "1"
            out *= current
        end
        current *= current
    end
    return out
end

function pow_mod(a::Polynomial, n::Int, p::Int)::Polynomial
    poly_mod = PolynomialModP(a,p)
    return pow_mod(poly_mod,n).poly
end
            



|