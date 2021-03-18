using PolymerArchitecture
using Test

using LightGraphs
using Polymer

include("chains.jl")

include("test_subtree.jl")
# include("test_equivalent.jl")
include("test_process.jl")

nothing # to avoid display last result of a @test block in REPL