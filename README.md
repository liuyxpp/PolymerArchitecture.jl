# PolymerArchitecture.jl

**PolymerArchitecture.jl** provides a *graph* representation of the polymer architecture. Based on such representation, equivalent sub-architectures and semi-equivalent sub-architectures of a non-cyclic block copolymer chain can be identified automatically. Equivalent sub-architectures and semi-equivalent sub-architectures are useful for reducing the computational cost of polymer field-theoretic simulations.

To construct a graph representation, polymers defined in the [Polymer.jl](https://github.com/liuyxpp/Polymer.jl) is required.

*Warning: Be aware that this package is currently under active development. The interface is highly unstable and subjects to change frequently.*

## Usage

### Installation

```julia
using Pkg
Pkg.add("PolymerArchitecture")
```

### Defining a Polymer

`PolymerArchitecture.jl` works with polymer models defined using [Polymer.jl](https://github.com/liuyxpp/Polymer.jl). First, define your polymer chain. Here is an example of an ABA triblock copolymer:

```julia
using Polymer
using Polymer: branchpoints, freeends

# Define segments
sA = KuhnSegment(:A)
sB = KuhnSegment(:B)

# Define topology points
eb = branchpoints(2)
fe = freeends(2)

# Define blocks
A1 = PolymerBlock(:A1, sA, 0.3, eb[1], fe[1])
A2 = PolymerBlock(:A2, sA, 0.3, eb[2], fe[2])
B = PolymerBlock(:B, sB, 0.4, eb[1], eb[2])

# Create the block copolymer
chainABA = BlockCopolymer(:ABA, [A1, A2, B])
```

### Building the Graph

Convert the `BlockCopolymer` object into a graph representation:

```julia
using PolymerArchitecture

g = BlockCopolymerGraph(chainABA)
```

### Classifying Architectures

You can identify the architecture type of the polymer chain:

```julia
arch_type = chaintype(chainABA)
println(arch_type) # Output: LinearArchitecture()
```

### Finding Equivalent Sub-Architectures

The package can identify equivalent blocks and sub-architectures, which is useful for symmetry reduction in polymer field-theoretic simulations (both SCFT and FTS).

```julia
eqblocks = group_equivalent_blocks(chainABA)
# Or
eqblocks = group_equivalent_blocks(g)
# 3-element Vector{Vector{Pair{Int64, Int64}}}:
#  [2 => 1, 4 => 3]
#  [1 => 3, 3 => 1]
#  [1 => 2, 3 => 4]
```

You can also consult the test codes reside in the `test` folder to learn more.

See [how PolymerArchitecture.jl is used in practice](https://www.yxliu.group/Polyorder.jl/dev/guide/guide1/#Create-Polymer-Systems).

## Contribute

* Star the package on [github.com](https://github.com/liuyxpp/PolymerArchitecture.jl).
* File an issue or make a pull request on [github.com](https://github.com/liuyxpp/PolymerArchitecture.jl).
* Contact the author via email <lyx@fudan.edu.cn>.

## Links

* [Yi-Xin Liu Research Group](http://www.yxliu.group)
* [Source code](https://github.com/liuyxpp/PolymerArchitecture.jl)
* [Additional resources](https://www.yxliu.group/Polyorder.jl/dev/guide/guide1/)