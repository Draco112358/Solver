using Test, GZip, JSON
include("../src/lib/utility.jl")


@testset let data = read_gzipped_json_file("/home/draco112358/esymiaProjects/mesherOutputs/init.gz")
    @test data isa Dict
    @test haskey(data, "mesher_matrices") == true
    @test haskey(data, "cell_size") == true
    @test haskey(data, "origin") == true
end