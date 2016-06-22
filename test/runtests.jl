using Hermetic
using DataStructures
using Base.Test



## lexicographic ordering for m = 3, k = 3

L =  [0     0     0
      0     0     1
      0     1     0
      1     0     0
      0     0     2
      0     1     1
      0     2     0
      1     0     1
      1     1     0
      2     0     0
      0     0     3
      0     1     2
      0     2     1
      0     3     0
      1     0     2
      1     1     1
      1     2     0
      2     0     1
      2     1     0
      3     0     0]


x = [0,0,0]
for j = 1:19
    mono_next_grlex!(x, 3)
    println("Checking glxr ordering: ", j+1)
    @test x == vec(L[j+1,:])
end

XX = zeros(Int, 20, 3);
mono_grlex!(XX, 3)
@test XX == L



He_table = OrderedDict([
              (1, [0, 1.]),
              (2, [-1., 0, 1]),
              (3, [0., -3, 0, 1]),
              (4, [3., 0, -6, 0, 1]),
              (5, [0., 15, 0, -10, 0, 1]),
              (6, [-15., 0, 45, 0, -15, 0, 1]),
              (7, [0, -105, 0, 105, 0, -21, 0, 1]),
              (8, [105, 0, -420, 0, 210, 0, -28, 0, 1]),
              (9, [0, 945, 0, -1260, 0, 378, 0, -36, 0, 1]),
              (10,[-945, 0, 4725, 0, -3150, 0, 630, 0, -45, 0, 1])
             ])

Hen_table = OrderedDict([
              (1, [0,1.]),
              (2, [-0.707107, 0, 0.707107]),
              (3, [0, -1.22474, 0, 0.408248]),
              (4, [0.612372, 0, -1.22474, 0, 0.204124]),
              (5, [0, 1.36931, 0, -0.912871, 0, 0.0912871]),
              (6, [-0.559017, 0, 1.67705, 0, -0.559017, 0, 0.0372678]),
              (7, [0, -1.47902, 0, 1.47902, 0, -0.295804, 0, 0.0140859]),
              (8, [0.522913, 0, -2.09165, 0, 1.04583, 0, -0.139443, 0, 0.00498012]),
              (9, [0, 1.56874, 0, -2.09165, 0, 0.627495, 0, -0.0597614, 0, 0.00166004]),
              (10,[-0.496078, 0, 2.48039, 0, -1.65359, 0, 0.330719, 0, -0.0236228, 0, 0.000524951])
                         ])

for i = enumerate(He_table)
    o, c, e = Hermetic.He_coefficients(i[1])
    @test i[2][:2][e + 1] == c
end

for i = enumerate(Hen_table)
    o, c, e = Hermetic.Hen_coefficients(i[1])
    @test_approx_eq_eps i[2][:2][e + 1] c 1e-05
end

## Hen_value


## Hen(j, 2) j = 0:10
Hen_2 = [1., 2., 2.1213203435596424, 0.8164965809277261,
         -1.0206207261596576, -1.6431676725154982, -0.4099457958749614, 
         1.2113877651108738, 1.2400496821844333, -0.31540754968546497,
         -1.3758956718849784]

## Hen(j, -2) j = 0:10
Hen_m2 = [1.,-2.,2.121320344,-0.8164965809,-1.020620726,1.643167673,
          -0.4099457959,-1.211387765,1.240049682,0.3154075497,-1.375895672]


## Hen(j, 1.2) j = 0:10
Hen_12 = [1.,1.2,0.3111269837,-0.7642407997,-0.7279883516,0.2928782059,
          0.8080398352,0.09533989072,-0.7154027646,-0.3760484168,0.5359903132]

## Hen(j, -1.2) j = 0:10
Hen_m12 = [1.,-1.2,0.3111269837,0.7642407997,-0.7279883516,-0.2928782059,
           0.8080398352,-0.09533989072,-0.7154027646,0.3760484168,0.5359903132]


for j = 0:10
    @test_approx_eq_eps [Hen_2[j+1]] Hermetic.Hen_value(j, [2.]) 1e-8
    @test_approx_eq_eps [Hen_m2[j+1]] Hermetic.Hen_value(j, [-2.]) 1e-8
    @test_approx_eq_eps [Hen_12[j+1]] Hermetic.Hen_value(j, [1.2]) 1e-8
    @test_approx_eq_eps [Hen_m12[j+1]] Hermetic.Hen_value(j, [-1.2]) 1e-8
end


## Polynomials add test

function poly_add_test()
    m  = 3
    o1 = 6
    c1 = [7.0, - 5.0, 9.0, 11.0, 0.0, - 13.0]
    e1 = [1, 2, 4, 5, 12, 33]

    o2 = 5
    c2 = [2.0, 3.0, -8.0, 4.0, 9.0]
    e2 = [1, 3, 4, 30, 33] 

    println("Add two polynomials")
    println("P(X) = P1(X) + P2(X)")
    println("")
    Hermetic.polynomial_print(m, o1, c1, e1, title = "P1(X) = ")
    Hermetic.polynomial_print(m, o2, c2, e2, title = "P2(X) = ")

    o, c, e = Hermetic.polynomial_add(o1, c1, e1, o2, c2, e2)
    Hermetic.polynomial_print(m, o, c, e, title = "P1(X) + P2(X) = ")
    @test o == 7
    @test c == [9.0, -5.0, 3.0, 1.0, 11.0, 4.0, -4.0]
    @test e == [1 2 3 4 5 30 33]

end

function poly_mul_test()
    m  = 3
    o1 = 4
    c1 = [2.0, 3.0, 4.0, 5.0]
    e1 = [  1,   3,   4,   6]

    o2 = 2
    c2 = [6.0, 7.0]
    e2 = [  2,   5]

    println("Multiply two polynomials")
    println("P(X) = P1(X) * P2(X)")
    println("")
    Hermetic.polynomial_print(m, o1, c1, e1, title = "P1(X) = ")
    Hermetic.polynomial_print(m, o2, c2, e2, title = "P2(X) = ")

    o, c, e = Hermetic.polynomial_mul(m, o1, c1, e1, o2, c2, e2)
    Hermetic.polynomial_print(m, o, c, e, title = "P1(X) * P2(X) = ")
    @test o == 7
    @test c == [8.0,9.0,9.0,10.0,21.0,11.0,12.0]
    @test e == [2, 5, 6, 8, 12, 15, 22]
end













