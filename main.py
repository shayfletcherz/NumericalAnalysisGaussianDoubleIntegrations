"""
Application for double integral calculation using the Gaussian quadrature methode.
Can be calculated polynomials up to order 9.
"""


import sympy as sym
from sympy import lambdify
import tableprint

x, y = sym.symbols('x y')

# Xi values (Legendre polynomials roots)

p = [[-0.5773502692, -0.7745966692, -0.8611363116, -0.9061798459],
     [0.5773502692, 0.0000000000, -0.3399810436, -0.5384693101],
     [0.0000000000, 0.7745966692, 0.3399810436, 0.0000000000],
     [0.0000000000, 0.0000000000, 0.8611363116, 0.5384693101],
     [0.0000000000, 0.0000000000, 0.0000000000, 0.9061798459]]

# Wi values
w = [[1.0000000000, 0.5555555556, 0.3478548451, 0.2369268850],
     [1.0000000000, 0.8888888889, 0.6521451549, 0.4786286705],
     [0.0000000000, 0.5555555556, 0.6521451549, 0.5688888889],
     [0.0000000000, 0.0000000000, 0.3478548451, 0.4786286705],
     [0.0000000000, 0.0000000000, 0.0000000000, 0.2369268850]]

"""
Calculate a particular integral using the Gaussian method in tow variable
Parameters: a,b - The boundaries of the integral, f- function, k - abscissas number
Return: F(y) function 
"""

def gaussian_quadrature_xy(a, b, f, k):

    sum = 0

    for i in range(5):

        t = (b - a) / 2 * p[i][k] + (b + a) / 2
        sum = sum + (f(t , y) * w[i][k])

    print("N:", k, ",F(y)=", (b - a) / 2 * sum)
    return (b - a) / 2 * sum

"""
Calculate a particular integral using the Gaussian method in one variable
Parameters: a,b - The boundaries of the integral, f- function, k - abscissas number
Return: Integral Approximate value
"""

def gaussian_quadrature_y(a, b, f, k):
    sum = 0
    for i in range(5):
        t = (b - a) / 2 * p[i][k] + (b + a) / 2
        sum = sum + (f(t) * w[i][k])

    return (b - a) / 2 * sum


"""
Calculate a particular integral using the Gaussian method in tow variable
Parameters: a,b - The boundaries of the integral, f- function, k - abscissas number
Return: Integral Approximate value
"""

def gaussian_quadrature_double_integral(f, a, b, c, d, k):

    fy = gaussian_quadrature_xy(a, b, f, k)
    fy = lambdify(y, fy)
    return (gaussian_quadrature_y(c, d, fy, k))

def main():

    tableprint.banner("Welcome to the Gaussian Quadrature Double Integral")

    #Choos the inner integral boundaries
    a = 0
    b = 2

    # Choos the outer integral boundaries
    c = -1
    d = 1

    #Choos abscissas number - Affects the level of accuracy
    k = 4

    #Set the polynomial
    fxy = x**7 + 2* x*y


    print("F(x,y)=", fxy, "\n")
    fxy = lambdify([x, y], fxy)

    data = []

    for i in range(k):

        data.append([i, gaussian_quadrature_double_integral(fxy, a, b, c, d, i)])

    headers = ['N number', 'Integral']

    tableprint.table(data, headers, '15g')


if __name__ == "__main__":
    main()

