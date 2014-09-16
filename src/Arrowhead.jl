module Arrowhead
using DoubleDouble

export GenSymArrow, SymArrow, aheig, bisect, invA, GenSymDPR1, SymDPR1, aheigall

include("arrowhead1.jl")
include("arrowhead3.jl")
end # module
