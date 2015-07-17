module Arrowhead
using DoubleDouble

export GenSymArrow, SymArrow, eig, bisect, inv, GenSymDPR1, SymDPR1, HalfArrow, GenHalfArrow, svd
# println("test")
include("arrowhead1.jl")
include("arrowhead3.jl")
include("arrowhead4.jl")
include("arrowhead5.jl")
include("arrowhead6.jl")

end # module
