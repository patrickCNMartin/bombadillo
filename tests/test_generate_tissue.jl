include("../src/MyModule.jl")
using .MyModule

@testset "Tissue tests" begin
    t = Tissue()
    @test t.cell_types === nothing
end
