# traits for the type of chain architecture
abstract type AbstractPolymerArchitecture end
# Non-cyclic architectures
abstract type NonCyclicArchitecture <: AbstractPolymerArchitecture end
struct LinearArchitecture <: NonCyclicArchitecture end
abstract type BranchedArchitecture <: NonCyclicArchitecture end
struct StarArchitecture <: BranchedArchitecture end
struct CombArchitecture <: BranchedArchitecture end
struct GeneralBranchedArchitecture <: BranchedArchitecture end
# Polymer that has ring(s) in it.
abstract type CyclicArchitecture <: AbstractPolymerArchitecture end
struct RingArchitecture <: CyclicArchitecture end

iscyclicchain(::CyclicArchitecture) = true
iscyclicchain(::AbstractPolymerArchitecture) = false
isnoncyclicchain(pa::AbstractPolymerArchitecture) = !iscyclicchain(pa)

"""
    islinearchain(<:AbstractPolymerArchitecture)

Check if the chain is linear. Currently we use explicit traits to check. It is possible to check the architecture by examining PolymerComponent instance. But it is not implemented yet.
"""
islinearchain(::LinearArchitecture) = true
islinearchain(::AbstractPolymerArchitecture) = false