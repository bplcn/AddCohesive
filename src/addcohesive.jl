module AddCohesive
    using AbaAccess
    using PBCHandler2D

    export addcohesive_2d_all!,addcohesive_2d!

    include("genernalfunctions.jl")
    include("addcohfunctions.jl")
    include("addcohfunctions3d.jl")



end