using NemoHyper
using Nemo
using Test

import Nemo: acb, AcbField, fmpq, QQ, libarb
import NemoHyper: nacb, compute_pfq

#### tests #####################################################################

function θ3(a::acb, b::acb)
  R = parent(a)
  z = [R(), R(), R(), R()]
  ccall((:acb_modular_theta, libarb), Nothing,
        (Ref{acb}, Ref{acb}, Ref{acb}, Ref{acb}, Ref{acb}, Ref{acb}, Int),
        z[1], z[2], z[3], z[4], a, b, precision(R))
  return z[3]
end

function gamma(a::Rational, C::AcbField)
  R = ArbField(precision(C))
  return C(gamma(QQ(a), R))
end

function hypergeometric_pfq(a::Vector, b::Vector, z::acb)
  C = parent(z)
  A = fmpq[QQ(i) for i in a]
  B = fmpq[QQ(i) for i in b]
  return C(compute_pfq(A, B, nacb(z), precision(C)))
end


function run_tests()
println()
println("***************************** tests *********************************")

#if @debug()
  @testset "initial values" begin
    CC = AcbField(100)
    for z in [CC(99//100), CC(101//100)], c in 1:5, k in 1:5
      @test isfinite(hypergeometric_pfq([1//3,2//3],[c],z))
      @test isfinite(hypergeometric_pfq([k+1//3,2//3],[c],z))
      @test isfinite(hypergeometric_pfq([1//3,k+2//3],[c],z))
    end
  end
#end

@testset "singular points" begin
  CC = AcbField(53)

  z = CC("[0.0 +/- 1.0e-10]")
  @test overlaps(hypergeometric_pfq([1,2,3],[4,5],z), one(CC))

  z = zero(CC)
  @test overlaps(hypergeometric_pfq([1,2,3],[4,5],z), one(CC))

  z = CC("[1.0 +/- 1.0e-10]")
  @test overlaps(hypergeometric_pfq([1,2,3],[4,5],z), 120-12*const_pi(CC)^2)

  z = one(CC)
  @test overlaps(hypergeometric_pfq([1,2,3],[4,5],z), 120-12*const_pi(CC)^2)

  z = CC("[1.0 +/- 1.0e-10]")
  @test !isfinite(hypergeometric_pfq([5,4,3],[2,1],z))

  z = one(CC)
  @test !isfinite(hypergeometric_pfq([5,4,3],[2,1],z))

end

@testset "algebraic cases" begin
  CC = AcbField(100)

  for j in [CC(1//7+2000), CC(1//7-2000),
            CC(1//7+1728), CC(1728), CC(1//7-1728),
            CC(1//7+200), CC(1//7-200),
            CC(1//7+64), CC(-64), CC(1//7-64),
            CC(1//7+2), CC(1//7-2)]
    f1 = hypergeometric_pfq([-1//12,1//4],[2//3],1728/j)
    f2 = hypergeometric_pfq([1//4, 7//12],[4//3],1728/j)
    t = j/3^3*f1^3/f2^3
    @test overlaps(j, 27*t*(t+8)^3/(t-1)^3)

    f1 = hypergeometric_pfq([-1//24,7//24],[3//4],1728/j)
    f2 = hypergeometric_pfq([5//24,13//24],[5//4],1728/j)
    t = j/2^4*f1^4/f2^4
    @test overlaps(j, 16*(t^2+14*t+1)^3/(t*(t-1)^4))

    f1 = hypergeometric_pfq([-1//60,19//60],[4//5],1728/j)
    f2 = hypergeometric_pfq([11//60,31//60],[6//5],1728/j)
    t = j*f1^5/f2^5
    @test overlaps(j, (t^4+228*t^3+494*t^2-228*t+1)^3/(t*(t^2-11*t-1)^5))

    f1 = hypergeometric_pfq([-1//42,13//42,9//14],[4//7,6//7],1728/j)
    f2 = hypergeometric_pfq([5//42,19//42,11//14],[5//7,8//7],1728/j)
    f3 = hypergeometric_pfq([17//42,31//42,15//14],[9//7,10//7],1728/j)
    t = j*f1*f2^2/f3^3
    @test overlaps(j, (t^2-t+1)^3*(t^6+229*t^5+270*t^4-1695*t^3+1430*t^2-
                                      235*t+1)^3/((t^2-t)*(t^3-8*t^2+5*t+1)^7))

    f1 = hypergeometric_pfq([-1//10,1//10,1//2],[1//5,4//5],-64/j)
    f2 = hypergeometric_pfq([1//10,3//10,7//10],[2//5,6//5],-64/j)
    f3 = hypergeometric_pfq([7//10,9//10,13//10],[8//5,9//5],-64/j)
    t = j*f2^2/f3*(j*f1^2 + f2*f3)/(j*f2^3 + f1*f3^2)
    @test overlaps(j, (t^2-1)*(t^2-4*t-1)^5/(t*(t^2+t-1)^5))
  end
end

@testset "function zeros" begin
  CC = AcbField(100)
  i = onei(CC)
  ipi = const_pi(CC)*i

  for m in [CC(1//10), CC(2//3), CC(9//10)]
    f1 = hypergeometric_pfq([1//3,2//3],[1],m)
    f2 = hypergeometric_pfq([1//3,2//3],[1],1-m)
    f3 = hypergeometric_pfq([5//6,6//6,7//6],[3//2,3//2],1-m)
    τ = sqrt(CC(-1//3))*f2/f1
    z = (2 + 3*τ)/4 - sqrt(1-m)*f3/(6*ipi*f1)
    @test contains_zero(θ3(z,2*τ)*θ3(z,6*τ) + exp(2*ipi*τ)*θ3(z+τ,2*τ)*θ3(z-3*τ,6*τ))
  end

  for j in [CC(10000//3), CC(10//3)]
    f1 = hypergeometric_pfq([1//12,5//12],[1],1728/j)
    f2 = hypergeometric_pfq([1//12,5//12],[1//2],1-1728/j)
    f3 = hypergeometric_pfq([1//3,2//3,1],[3//4,5//4],1-1728/j)
    τ = i*(2*gamma(1//2,CC)/gamma(7//12,CC)/gamma(11//12,CC)*f2/f1 - 1)
    z = (1 + τ)/2 + sqrt(CC(2//3))/ipi*(1-1728/j)^(1//4)*f3/f1
    @test contains_zero(ellipwp(z, τ))
  end
end

@testset "entire case" begin
  CC = AcbField(100)

  for z in [CC(1), CC(-1), CC(10), CC(-10)]
    @test overlaps(hypergeometric_pfq([1],[2],z), (exp(z)-1)/z)
  end
end

@testset "divergent case" begin
  CC = AcbField(100)

  # 0! - 1! + 2! - 3! + ... = 0.59...
  @test overlaps(hypergeometric_pfq([1,1],[],CC(-1)),
                 CC("[0.59634736232319407434107849937 +/- 2.01e-30]"))

  # 0!^2 - 1!^2 + 2!^2 - 3!^2 + ... = 0.66...
  @test overlaps(hypergeometric_pfq([1,1,1],[],CC(-1)),
                 CC("[0.66809132637777776543301498750 +/- 1.58e-30]"))

  # 0!^3 - 1!^3 + 2!^3 - 3!^3 + ... = 0.72...
  @test overlaps(hypergeometric_pfq([1,1,1,1],[],CC(-1)),
                 CC("[0.72374994418155305222032415740 +/- 5.33e-30]"))
end

end

function run_bench(a, b, z, p, fx = "", fy = "")

  A = map(QQ, a)
  B = map(QQ, b)
  println("---- $A $B $z digits=$p ----")

  p = trunc(Int, p*3.322)
  zz = nacb(narb(QQ(real(z)), p), narb(QQ(imag(z)), p))
  @time f = compute_pfq(A, B, zz, p)
  @assert isempty(fx) || overlaps(set!(narb(), fx, p), real(f))
  @assert isempty(fy) || overlaps(set!(narb(), fy, p), imag(f))
  bits_lost = p - precision(f)
  return bits_lost
end

function run_benchmarks()
println()
println("************************** benchmarks *******************************")

  t1 = time()
  bits_lost = 0

  bits_lost += run_bench([1//2,1//3,1//5,1//4],[8//7,1//6,1//9], 1//10, 20, "[1.04165510107809208771 +/- 2.43e-21]", "[+/- 2.07e-24]")
  bits_lost += run_bench([1//2,1//3,1//5,1//4],[8//7,1//6,1//9],-1+im//300, 200, "[0.7361464107568699653654181932021195729016525520251669994612965514688945567873978145804140256340981534032171785111170341312718204804372620880207341602115123749944939921894542469460886808720225689828743232125 +/- 7.05e-206]", "[0.0006091142434506593619507708486086804762150980760156332916347162351588534898343587334997033663903136796976875663495423707271793534625090625898822976732269967397826684043934917897380440587935220808422516 +/- 1.57e-203]")
  bits_lost += run_bench([1//2,1//3,1//5,1//4],[8//7,1//6,1//9],1//300+im, 200, "[0.853498211060579954490363364784481268197877903709746131856530610099120791429240016786586006359060701719614321941609976993732551828969874878632105881859561180925930687147324142274043333887472085663429960 +/- 4.43e-202]", "[0.302235718305177081880438308463901805943967344220865402810895162185350037204941683860825775142412188198314838346748247659487319633322537866130218549109342613855275799727201394499484865564020106350734127 +/- 2.99e-202]")
  bits_lost += run_bench([1//2,1//3,1//5,1//4],[8//7,1//6,1//9], 11//10, 100, "[2.09473837006811761028534487482955095697725399923738208299473100086979606619943028860894705913795820 +/- 2.23e-99]", "[-1.32063424488884879145689203640820495365830310353133670921975318942279074195339938129952186569083161 +/- 3.11e-99]")
  bits_lost += run_bench([1//2,1//3,1//5,1//4],[8//7,1//6,409//420], 11//10, 100, "[1.10012992690145570229907705218466524693847390580609878681247630305435529779555060922581931979395195 +/- 7.02e-100]", "[-0.019650631936840220248345681395791981619880186753148199866219643526635725512689324474868785073820769 +/- 1.99e-100]")
  bits_lost += run_bench([1//2,1//3,1//5,1//4],[8//7,1//6,1//9], 21//10+im, 100, "[0.71474921490254889156840981949000811062050067188575315042564339656971171954121748902674586037861411 +/- 4.77e-99]", "[0.9781249447537276082866697472831284801654120848168295066054521588403569889034210483850430568298298 +/- 3.56e-98]")
  bits_lost += run_bench([1//2,1//3,1//5,1//4],[8//7,1//6,10//9], 21//10, 100, "[1.09611960714180653165793929034538750471702782740589204327790995048489804128469673421588286818944088 +/- 2.55e-99]", "[-0.1103381390329406104849629727151339938760273060221215999605687131539917991045753083253694103910355 +/- 4.15e-98]")
  bits_lost += run_bench([1//2,1//3,3//2,3//2],[8//7,1//6,409//420], 21//10*im, 100, "[-0.3836936156357715109532470008774342030087022878339141902193492463878637848072764973662114392423830 +/- 3.82e-98]", "[0.3510560688879407205546122233658260321782263643876690384980560383966660407148203677984173909952239 +/- 4.45e-98]")
  bits_lost += run_bench([1//2,1//3,3//2,3//2],[8//7,1//6,7//2+1//4], 21//10*im, 100, "[0.498286277362499154454672287676933977606853722206194319858084629633516142945046940166960257970613 +/- 3.22e-97]", "[0.5388483851104630403622090287525048639997504587807914458817597010957463072020568473516818023683285 +/- 4.40e-98]")
  bits_lost += run_bench([1//2,1//3,3//2,3//2],[8//7,1//6,409//420], 21//10, 100, "[-0.3778888881782393774553127309816076760288807629794011361623259687233181125952890271460921393790923 +/- 1.00e-98]", "[1.0324057197731830866862412722661912510450592489443071095814531303361891341977524868140295288542887 +/- 3.18e-98]")
  bits_lost += run_bench([1//2,1//3,3//2,3//2],[8//7,1//6,7//2+1//4], 21//10, 100, "[0.5587476752887922757498928650974538228321146670224016380361419696854589138256157940556936141300072 +/- 7.32e-98]", "[-2.649787031624682836178763161718415233309248545654569640464561282640046441505973680884328801960931 +/- 1.55e-97]")
  bits_lost += run_bench([1//2,1//3,3//2,5//2],[8//7,1//6,7//2+1//4], 21//10, 100, "[-1.799478585201052861552266948753576293798081804909492183540383822249557260558139389282027739306081 +/- 2.29e-97]", "[-2.0550213095299542518482540717224091776461005489751204129826789474740506081260337268088507033546378 +/- 5.82e-98]")
  bits_lost += run_bench([1//2,1//3,1//3,3//2,3//2],[8//7,1//6,7//2+1//5,1//4], 21//10, 100, "[-0.116857050064488324880005256242593338635544257506370503116752352065017363357264920954229325317036 +/- 5.51e-97]", "[-3.515406819990606948896449062181976165047281121147375761111710596644569438443262697365115597303044 +/- 2.31e-97]")
  bits_lost += run_bench([1//2,1//3,1//3,3//2,3//2,4//3],[8//7,1//6,7//2,1//4,1//5], 21//10, 100, "[-19.9404485154538633201163293612929060982693170333095030435564215171671004250390889987049810825485 +/- 2.54e-95]", "[17.0022670479727611109224091973883810218739134352957349354752346919984236719406436188027263944349 +/- 4.04e-95]")
  bits_lost += run_bench([1//2,1//3,1//3,3//2,3//2,4//3],[8//7,1//6,7//2,1//4,5+1//5], 101//100, 100)
  bits_lost += run_bench([1,1,2,2,3],[3//2,3//2,3//2,3//2], 21//10, 100, "[-0.51585621412817788953636435138740344792545648221588737349075642969025142618998460938638704600278106 +/- 3.91e-99]", "[-0.05470891931824783630344691907005460888039986744824234445189809615982014632842245019368879129590597 +/- 3.12e-99]")
  bits_lost += run_bench([-2,1,2,2,3],[4//3,4//3,4//3,4//3], 21//10, 100, "[25.724630102040816326530612244897959183673469387755102040816326530612244897959183673469387755102041 +/- 2.05e-97]")

  println("total: ", time()-t1)
  println("bits lost: ", bits_lost)
end

function run_heuristics()
println()
println("************************** heuristics *******************************")
  #if @debug()
    #fac = 1
  #else
    fac = 100
  #end
  for (z, limit) in [(1, 250*fac), (1-1//50, 60*fac)]
    d = 20
    while (d = 2*d) <= limit
      m = trunc(Int, d*4.0)
      N = trunc(Int, d*1.0)
      p = trunc(Int, d*3.322)
      println()
      println("--------- $d digits, z = $z ----------")

      @time s1 = heuristic_hyp(map(QQ,[1//2,1//3,1//4,3//5]),
                               map(QQ, [11//7,7//6,1//8]), QQ(z),
                               m, N, p)

      @time s2 = hypergeometric_pfq([1//2,1//3,1//4,3//5],
                                    [11//7,7//6,1//8], AcbField(p)(z))
      s2 = nacb(s2)
      @show sub(s1, s2, p+60)
      if !overlaps(s1, s2)
        println("oops: ", sub(s2, s1, p+100))
        error("bad answer")
      end
    end
  end
end

run_tests()
run_benchmarks()
run_benchmarks()
run_heuristics()



println()
println("************************ bad borel sum ****************************")

function rescale(Aθ::Vector, s)

  Fθ = parent(Aθ[1])
  θ = gen(Fθ)
  F = base_ring(Fθ)
  Fx, x = PolynomialRing(F, "x")
  Ax = transpose_vars(Aθ, Fx)

  x = gen(Fx)
  ra = length(Ax) - 1
  eq = elem_type(Fx)[Ax[1+ra]]
  for k in ra-1:-1:0
    eq = elem_type(Fx)[
          ((i<length(eq)) ? x*derivative(eq[1+i]) : zero(Fx)) +
          ((i > 0) ? eq[1+i-1] : Ax[1+k])
        for i in 0:length(eq)]
  end

  eq = map(p -> evaluate(p, s*x), eq)

  # cancel content while θ is on the right
  c = map(p -> p, eq)
  cont = zero(Fx)
  for ci in c
    cont = gcd(cont, ci)
  end
  c = elem_type(Fx)[divexact(ci, cont) for ci in c]

  # move θ to the left
  Cx = elem_type(Fx)[]
  for ci in reverse(c)
    equ_mul_θ_add!(Cx, ci, Fx)
  end
  return transpose_vars(Cx, Fθ)
end


# a0(x)*(xf) + a1(x)*θ(xf) + a2(x)*θθ(xf) = 0
# θg = xf
# θθg = θxf
# int(f) = g
function integrate(Aθ::Vector)
  Fθ = parent(Aθ[1])
  θ = gen(Fθ)
  F = base_ring(Fθ)
  Fx, x = PolynomialRing(F, "x")
  Ax = transpose_vars(Aθ, Fx)

  Fx = parent(Ax[1])
  x = gen(Fx)
  ra = length(Ax) - 1
  eq = elem_type(Fx)[Ax[1+ra]]
  for k in ra-1:-1:0
    eq = elem_type(Fx)[
          ((i<length(eq)) ? x*derivative(eq[1+i])-eq[1+i] : zero(Fx)) +
          ((i > 0) ? eq[1+i-1] : Ax[1+k])
        for i in 0:length(eq)]
  end

  pushfirst!(eq, zero(Fx))

  # cancel content while θ is on the right
  c = map(p -> p, eq)
  cont = zero(Fx)
  for ci in c
    cont = gcd(cont, ci)
  end
  c = elem_type(Fx)[divexact(ci, cont) for ci in c]

  # move θ to the left
  Cx = elem_type(Fx)[]
  for ci in reverse(c)
    equ_mul_θ_add!(Cx, ci, Fx)
  end
  return transpose_vars(Cx, Fθ), Cx
end

function shift(Aθ::Vector, s)
  Fθ = parent(Aθ[1])
  θ = gen(Fθ)
  F = base_ring(Fθ)
  Fx, x = PolynomialRing(F, "x")
  Ax = transpose_vars(Aθ, Fθ)

  FFx = FractionField(Fx)

  eq = elem_type(FFx)[]
  for k in length(Aθ):-1:1
    # y = x + s
    # θ_x = (1-s/y) θ_y
    # apply θ_x and add Ax[k]
    # (1-s/y) θ (Σ_i c_i(y)*θ^i(g))
    # = (1-s/y) θ (Σ_i c_i(y)*θ^i(g))
    # = Σ_i (1-s/y) y d/dy (c_i(y)*θ^i(g))
    # = Σ_i (1-s/y) (y*c_i'(y)*θ^i(g) + c_i(y)*θ^{i+1}(g))
    eq = elem_type(FFx)[(
          ((i<length(eq)) ? (x+s)*derivative(eq[1+i]) : zero(FFx)) +
          ((i > 0) ? (x+s)//x*eq[1+i-1] : FFx(evaluate(Ax[k], x+s))))
        for i in 0:length(eq)]
  end

  # cancel content while θ is on the right
  c = map(p -> (numerator(p), denominator(p)), eq)
  cont = zero(Fx)
  den = one(Fx)
  for ci in c
    cont = gcd(cont, ci[1])
    den = lcm(den, ci[2])
  end
  c = elem_type(Fx)[divexact(ci[1], cont)*divexact(den, ci[2]) for ci in c]

  # move θ to the left
  Cx = elem_type(Fx)[]
  for ci in reverse(c)
    equ_mul_θ_add!(Cx, ci, Fx)
  end
  return transpose_vars(Cx, Fθ), Cx
end



let T = fmpq, Φ = 200, CC = AcbField(Φ), a = T[QQ(1), QQ(1)], b = T[]
  zz = fmpq(-1//10)
  z = CC(zz)

println("a = $a   b = $b   z = $zz")


  F = parent(a[1])
  Fx, x = PolynomialRing(F, "x")
  # 0 -> z
  Pθ, Px = hyp_equ(a, b, Fx(0)//Fx(1), x//Fx(1))
  Fθ = parent(Pθ[1])
  θ = gen(Fθ)

  # eq for f(t)
  eq = hadamard_product(Pθ, [θ, Fθ(-1)])
  eq = θ .* eq  # hmmm.... bad corner case

  # eq for f(z*t)
  eq = rescale(eq, zz)

  # eq for exp(-t)*f(z*t)
  eq = cauchy_product(eq, [θ, Fθ(1)])

  # eq for int exp(-t)*f(z*t) dt
  # singular points at t = 0, t = 1/z, and t = inf
  Pθ, Px = integrate(eq)

  zi = 2
  e = eval_basis(Pθ, Px, zBranchInfo{T}(F(0), F(2)), [fmpz(0), fmpz(1)], [F(0), F(1)], 2, Φ)
println("integral from 0 to $zi:")
@show e[1,2]
@show e[2,2]

  while zi < 150
    nzi = zi + 1 + cld(zi, 3)
    nPθ, nPx = shift(Pθ, F(zi))
    e1 = eval_basis(nPθ, nPx, zBranchInfo{T}(F(0), F(nzi-zi)), [fmpz(0), fmpz(1)], [F(0), F(1)], 2, Φ)
    e = mul(e1, e, Φ+20)
    zi = nzi
println("integral from 0 to $zi:")
@show e[1,2]
@show e[2,2]
  end

println("real answer:")
@show hypergeometric_pfq(a,b,z)
