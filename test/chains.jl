branchpoints(n, prefix="EB") = [BranchPoint(Symbol(prefix*string(i))) for i in 1:n]

freeends(n, prefix="A") = [FreeEnd(Symbol(prefix*string(i))) for i in 1:n]

function chainAB3A3()
    sA = KuhnSegment(:A)
    sB = KuhnSegment(:B)
    eb0 = BranchPoint(:EB0)
    eb1 = BranchPoint(:EB1)
    eb2 = BranchPoint(:EB2)
    eb3 = BranchPoint(:EB3)
    A = PolymerBlock(:A, sA, 0.1, FreeEnd(:A), eb0)
    B1 = PolymerBlock(:B1, sB, 0.2, eb1, eb0)
    B2 = PolymerBlock(:B2, sB, 0.2, eb2, eb0)
    B3 = PolymerBlock(:B3, sB, 0.2, eb3, eb0)
    A1 = PolymerBlock(:A1, sA, 0.1, eb1, FreeEnd(:A1))
    A2 = PolymerBlock(:A2, sA, 0.1, eb2, FreeEnd(:A2))
    A3 = PolymerBlock(:A3, sA, 0.1, eb3, FreeEnd(:A3))
    return BlockCopolymer(:AB3A3, [A, B1, B2, B3, A1, A2, A3])
end

function branchAB()
	sA = KuhnSegment(:A)
	sB = KuhnSegment(:B)
	eb = branchpoints(8)
	fe = freeends(6, "B")
	A1 = PolymerBlock(:A1, sA, 0.12, eb[1], eb[3])
	A2 = PolymerBlock(:A2, sA, 0.12, eb[2], eb[3])
	A3 = PolymerBlock(:A1, sA, 0.12, eb[6], eb[8])
	A4 = PolymerBlock(:A1, sA, 0.12, eb[6], eb[7])
	A5 = PolymerBlock(:A1, sA, 0.2, eb[4], eb[5])
	B1 = PolymerBlock(:B1, sB, 0.02, eb[1], fe[1])
	B2 = PolymerBlock(:B2, sB, 0.02, eb[1], fe[2])
	B3 = PolymerBlock(:B3, sB, 0.02, eb[2], fe[3])
	B4 = PolymerBlock(:B4, sB, 0.02, eb[8], fe[4])
	B5 = PolymerBlock(:B5, sB, 0.02, eb[7], fe[5])
	B6 = PolymerBlock(:B6, sB, 0.02, eb[7], fe[6])
	B7 = PolymerBlock(:B7, sB, 0.08, eb[3], eb[4])
	B8 = PolymerBlock(:B8, sB, 0.12, eb[5], eb[6])
	return BlockCopolymer(:AB3A3, [A1, A2, A3, A4, A5, B1, B2, B3, B4, B5, B6, B7, B8])
end

abg = diblock_chain() |> BlockCopolymerGraph

starAB3A3g = chainAB3A3() |> BlockCopolymerGraph

branchg = branchAB() |> BlockCopolymerGraph