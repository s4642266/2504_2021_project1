
include("poly_factorization_project.jl")
println("For this repository, a new type was made for polynomials.")
println("The following are tests for the different functions that can be done with this new type")


x = x_poly()
@show p1 = 2x^3 + 4x^2 - 3x
@show p2 = 2x^4 - 4x^2 - 3x + 3

println("The first test is for polynomial mutliplication")
@show p1*p2
println("The following is a test for polynomial power")
@show p1^3

println("In the following 2 lines, it can be seen that the derivative function works on the polynomial class")
println("The first line is the derivative function on the polynomials multiplied togeher")
println("The second line is using the product rule")
d1 = derivative(p1*p2)
d2 = derivative(p1)*p2 + p1*derivative(p2);
println("Derivative of multiplied function: ", d1)
println("Derivative with product rule: ", d2)

prime = 17
p = mod((7x^3 + 2x^2 + 8x + 1)*(x^2+x+1),prime)
println("Will factor this polynomial (mod $prime): ", p)
factorization = factor(p,prime)
println("Here is the factorization: ", factorization)

pr = mod(expand_factorization(factorization),prime)
println("Reconstructing: ", pr)

println("The following is a test for polynomial addition")
x = x_poly()
@show p1 = 2x^2+3x +(-5) #Note we need +(-5)... try later without... how to fix?
@show p2 = -3x^2 - 4x +6
@show p1+p2

extended_euclid_alg(p1*p2,p2,101)