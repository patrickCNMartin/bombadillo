# this part include the other files you need 

using ArgParse
using FileIO


# example of a command line parser 
# taken from https://github.com/roberthoenig/WaveFunctionCollapse.jl/tree/master
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "coordinates"
            required = true
        "counts"
            required = true
        "--iter"
            arg_type = Int
            default = 10000
        
    end

    return parse_args(s, as_symbols=true)
end

function main()
    parsed_args = parse_commandline()
    filename = parsed_args[:filename]
    delete!(parsed_args, :filename)
    img = generate(filename; parsed_args...)
    FileIO.save("output.png", img)
    println("Done! Saved the output image as output.png.")
end

main()