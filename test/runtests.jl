using Barycentric
using LinearAlgebra
using Test

# These tests arbitrarily use 5 and 6 degree polynomials (odd and even are used to
# ensure coverage)

equi₅ = Equispaced{5}()
equi₆ = Equispaced{6}()
cheb1₅ = Chebyshev1{5}()
cheb1₆ = Chebyshev1{6}()
cheb2₅ = Chebyshev2{5}()
cheb2₆ = Chebyshev2{6}()

x5 = [0.814723686393179,1.7205156234688,1.8475024397623,2.76087829590132,3.39323754212673,3.49077794712614]
x6 = [0.278498218867048,0.825379738072032,1.78288657350633,2.74777510870561,2.90538819038315,3.87598097214377,4.83314792038672]
arb₅ = ArbitraryPolynomial(x5)
arb₆ = ArbitraryPolynomial(x6)

@testset "Nodes" begin
    @test nodes(equi₅) ≈ [-1,-0.6,-0.2,0.2,0.6,1]
    @test nodes(equi₆) ≈ [-1,-0.666666666666667,-0.333333333333333,0,0.333333333333333,0.666666666666667,1]

    @test nodes(cheb1₅) ≈ -[0.965925826289068,0.707106781186548,0.258819045102521,-0.258819045102521,-0.707106781186547,-0.965925826289068]
    @test nodes(cheb1₆) ≈ -[0.974927912181824,0.78183148246803,0.433883739117558,6.12323399573677e-17,-0.433883739117558,-0.781831482468029,-0.974927912181824]

    @test nodes(cheb2₅) ≈ -[1,0.809016994374947,0.309016994374947,-0.309016994374947,-0.809016994374947,-1]
    @test nodes(cheb2₆) ≈ -[1,0.866025403784439,0.5,6.12323399573677e-17,-0.5,-0.866025403784439,-1]

    @test nodes(arb₅) ≈ x5
    @test nodes(arb₆) ≈ x6
end

@testset "Weights" begin
    @test weights(equi₅) ≈ [-1,5,-10,10,-5,1]
    @test weights(equi₆) ≈ [1,-6,15,-20,15,-6,1]

    @test weights(cheb1₅) ≈ -[0.258819045102521,-0.707106781186547,0.965925826289068,-0.965925826289068,0.707106781186548,-0.258819045102521]
    @test weights(cheb1₆) ≈ -[0.222520933956314,-0.623489801858733,0.900968867902419,-1,0.900968867902419,-0.623489801858734,0.222520933956314]

    @test weights(cheb2₅) ≈ -[0.5,-1,1,-1,1,-0.5]
    @test weights(cheb2₆) ≈ -[0.5,-1,1,-1,1,-1,0.5]

    @test weights(arb₅) ≈ [-0.079601734468871,2.82206348049722,-3.28654127480991,1.17155138701322,-2.4317744194857,1.80430256125404]
    @test weights(arb₆) ≈ [0.0114361913589559,-0.0390633084690722,0.100394094241955,-0.588770483980979,0.552871013562898,-0.0415349717143145,0.00466746500055649]
end

xnew = [-0.514624351277159,0.285656117611641,0.427542456238857,0.849303738865132] # randomly generated

@testset "Interpolation matrix" begin
    @test interpolation_matrix(cheb1₅, xnew) ≈ [-0.0289845478236897 0.0959627734159381 -0.207065941756692 0.626076869404405 0.609098232472457 -0.0950873857124182;-0.0105623302174427 0.0465782216908936 0.999201970730136 -0.0492504664603377 0.0197735254445248 -0.00574092118777367;-0.0706583814898062 0.371759729030326 0.841448464571141 -0.206847346420822 0.0915972568975815 -0.0272997225884207;0.362858671520288 0.813049119299893 -0.267459002051222 0.142520710886204 -0.0742818875343108 0.0233123878791484][:, end:-1:1]
    @test interpolation_matrix(cheb1₆, xnew) ≈ [0.0127857533810382 -0.0411607067288632 0.0812980146502686 -0.166310956920905 0.955056229650331 0.199706640020739 -0.0413749740526084;0.041385687119839 -0.161088285741265 0.779202698741526 0.448772538954965 -0.160518366877509 0.0748749109083769 -0.0226291831059332;0.00285521790063521 -0.0123604217994083 0.997916710684666 0.0164279128299457 -0.00734604100263827 0.00362101552128418 -0.00111439413448407;0.184921983997812 0.964706061561875 -0.226419209876787 0.122921499942711 -0.0733011112746609 0.0399052720352278 -0.0127344963861772][:, end:-1:1]
    @test interpolation_matrix(cheb1₅, nodes(cheb1₅)) ≈ I
    @test interpolation_matrix(cheb1₆, nodes(cheb1₆)) ≈ I
end

@testset "Interpolate" begin
    c1 = Chebyshev1{10}()
    @test interpolate(c1, sin.(nodes(c1)), xnew) ≈ sin.(xnew)
    c2 = Chebyshev2{20}(0, 10)
    x = 10*rand(5)
    @test interpolate(c2, sin.(nodes(c2)), x) ≈ sin.(x)
end

@testset "Differentiation matrix" begin
    @test differentiation_matrix(cheb1₅) ≈ -[7.20976852010751 -10.5558337350587 5.27791686752937 -3.04720672422855 1.63299316185545 -0.517638090205042;1.4142135623731 0.707106781186547 -3.04720672422855 1.4142135623731 -0.707106781186548 0.218779599482357;-0.378937381963012 1.63299316185545 0.138700708242031 -1.93185165257814 0.757874763926024 -0.218779599482357;0.218779599482357 -0.757874763926024 1.93185165257814 -0.138700708242031 -1.63299316185545 0.378937381963012;-0.218779599482357 0.707106781186547 -1.4142135623731 3.04720672422855 -0.707106781186545 -1.4142135623731;0.517638090205041 -1.63299316185545 3.04720672422854 -5.27791686752936 10.5558337350587 -7.2097685201075]
    @test differentiation_matrix(cheb1₆) ≈ -[9.84466088119818 -14.5105621059791 7.48352452737972 -4.60952974192497 2.87399478545474 -1.59494677776481 0.512858431636276;1.84827792218219 1.0055981139743 -4.15304279314467 2.05143372654511 -1.18863516903897 0.639524003844967 -0.203155804362921;-0.4564869487803 1.98886610088761 0.267253750914115 -2.55809601537987 1.15238243548124 -0.569229914492255 0.17531059136945;0.22824347439015 -0.797473388882404 2.07652139657234 -1.02695629777827e-15 -2.07652139657234 0.797473388882405 -0.22824347439015;-0.175310591369451 0.569229914492254 -1.15238243548124 2.55809601537987 -0.26725375091411 -1.98886610088761 0.456486948780299;0.203155804362921 -0.639524003844966 1.18863516903897 -2.05143372654511 4.15304279314467 -1.00559811397431 -1.84827792218218;-0.512858431636278 1.59494677776481 -2.87399478545474 4.60952974192498 -7.48352452737972 14.5105621059791 -9.84466088119817]
    c1 = Chebyshev2{20}(0, 10)
    x = 10*rand(5)
    D = differentiation_matrix(c1)
    @test interpolate(c1, D*sin.(nodes(c1)), x) ≈ cos.(x)
end
