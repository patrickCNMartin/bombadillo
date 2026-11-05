@testset "network construction" begin
    edges = [Edge("a", "b"; sign = 1, K = 10),
             Edge("c", "b"; sign = -1, module_name = "m2"),
             Edge("L", "a"; sign = 1, layer = :chromatin, module_name = "m2")]
    net = RegNetwork(edges; ligands = ["L"])
    @test net.genes == ["a", "b", "c"]
    @test n_genes(net) == 3 && B.n_ligands(net) == 1
    bnd, chr = net.binding, net.chromatin
    # CSC: gene b (2) has both binding edges, sorted by source
    @test bnd.colptr == Int32[1, 1, 3, 3]
    @test bnd.src == Int32[1, 3] && bnd.sign == Int8[1, -1]
    # ligand inputs are indexed after genes
    @test chr.colptr == Int32[1, 2, 2, 2] && chr.src == Int32[4]
    @test net.module_genes == [[1, 2], [1, 2, 3]]

    progs = Programs(net, ["only_core" => ["core"], "both" => ["core", "m2"]])
    @test progs.mask == UInt64[0b01, 0b11]
    @test progs.A0[:, 1] == Float32[1, 1, 0]
    @test progs.A0[:, 2] == Float32[1, 1, 1]

    @test_throws ArgumentError Edge("a", "b"; sign = 0)
    @test_throws ArgumentError Edge("a", "b"; layer = :rna)
    @test_throws ArgumentError RegNetwork([Edge("a", "b"), Edge("a", "b")])
    @test_throws ArgumentError RegNetwork([Edge("a", "L")]; ligands = ["L"])
    @test_throws ArgumentError RegNetwork([Edge("x", "b")]; genes = ["b"])
    @test_throws ArgumentError Programs(net, ["p" => ["nope"]])
    @test_throws ArgumentError Kinetics(3; k_off = 1.5)
    @test_throws ArgumentError repressilator_network(genes = ["a", "b", "c", "d"])
end
