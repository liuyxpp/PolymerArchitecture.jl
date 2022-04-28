using PolymerArchitecture
using Test

using Graphs
using Polymer

include("chains.jl")

include("test_graph.jl")
include("test_subtree.jl")
include("test_equivalent.jl")
include("test_process.jl")
include("test_semi_equivalent.jl")

nothing # to avoid display last result of a @test block in REPL