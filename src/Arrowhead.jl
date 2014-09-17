module Arrowhead
using DoubleDouble

export GenSymArrow, SymArrow, aheig, bisect, invA, GenSymDPR1, SymDPR1, aheigall

include("arrowhead1.jl")
include("arrowhead3.jl")
include("arrowhead4.jl")

end # module
