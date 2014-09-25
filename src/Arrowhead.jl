module Arrowhead
using DoubleDouble

export GenSymArrow, SymArrow, aheig, bisect, invA, GenSymDPR1, SymDPR1, aheigall, dpr1eig, dpr1eigall, HalfArrow, GenHalfArrow, ahsvd, ahsvdall
# println("test")
include("arrowhead1.jl")
include("arrowhead3.jl")
include("arrowhead4.jl")
include("arrowhead5.jl")
include("arrowhead6.jl")

end # module
