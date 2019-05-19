module Arrowhead

    __precompile__(true)
    # using DoubleDouble
    using LinearAlgebra
    using Polynomials

    export GenSymArrow, SymArrow, eigen, bisect, inv, GenSymDPR1, SymDPR1, HalfArrow, GenHalfArrow, svd, tdc, rootsah, rootsWDK

    include("DoubleDouble.jl")
    include("arrowhead1.jl")
    include("arrowhead3.jl")
    include("arrowhead4.jl")
    include("arrowhead5.jl")
    include("arrowhead6.jl")
    include("arrowhead7.jl")
    include("arrowhead8.jl")
end # module
