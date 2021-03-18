using PolymerArchitecture
using Test

using LightGraphs
using Polymer

include("chains.jl")

include("test_subtree.jl")
include("test_equivalent.jl")

nothing # to avoid display last result of a @test block in REPL