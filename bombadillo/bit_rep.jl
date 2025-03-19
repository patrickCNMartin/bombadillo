using Random
using CSV
using DataFrames

# Let's get the data from the data directory


function bit_rep(vec)
    vec = join(string.(vec))
    bit_array = parse(BigInt,vec, base = 2)
    return bit_array
end

function array_to_bvec(array)
    vec = parse(BitVector,array, base = 2)
    return vec
end

function vec_rep(bit)
    vec_array = digits(bit, base = 2)
    return vec_array
end

function binarize_genes(dat, threshold = 0.7)
   dat[dat > quantile(dat, threshold),:] = 1
   dat[dat <= quantile(dat, threshold),:] = 0
   return dat
end

function binarize_matrix(data, threshold = 0.7)
    for row in Tables.namedtupleiterator(data)
        data[row,2:end] = binarize_genes(data[row,2:end], threshold)
    end
    return data
end

function binarize_neighbors(counts, graph)

end

# data_dir = "/Users/martinp4/Documents/exercise_data/"
# files = readdir(data_dir, join=true)
# coord = CSV.File(files[4]) |> DataFrame
# counts = CSV.File(files[5]) |> DataFrame


