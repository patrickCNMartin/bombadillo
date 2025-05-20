

function get_ecosystem(
    local_niche::Vector{Int},
    dist::Vector{Int};
    cell_labels::Vector{Int64},
    types::Vector,
    wave_diffusion::Int64)::AbstractMatrix
    ecosystem =  zeros(length(types), wave_diffusion)
    for cell in eachindex(local_niche)
        depth_index = dist[cell]
        depth_index == 0 && continue
        cell_index = cell_labels[local_niche[cell]]
        ecosystem[cell_index, depth_index] += 1
    end
    return ecosystem
end

function aggregate_shift(
    shift_set::Int64,
    ecosystem::Vector{Float64},
    shift_type::String,
    density_damp::Float64,
    wave_damp::Float64)
    if shift_type == "type_shifts"
        shift = shift_set
    elseif shift_type == "contact_shifts"
        density = ecosystem[1]
        shift = density_decay(density, shift_set, density_damp)
    elseif shift_type == "wave_shifts"
        shift = wave_decay(ecosystem, shift_set, wave_damp, density_damp)
    else
        throw(DomainError("Unrecognized type shift - Select from type_shifts, contact_shifts, wave_shifts"))
    end
    return shift
end




function apply_shift(
    gene_vec::Vector{Int64},
    gene::Int64,
    shift::Int64) 
    if gene + shift < 1
        shift = -gene + 1
    end
    if gene + shift > maximum(gene_vec)
        shift =  maximum(gene_vec) - gene
    end
    gene_vec[gene] = gene + shift
    return gene_vec
end



function ecosystem_shifts(
    base_cells::BaseCells,
    cell_label::Int64,
    ecosystem::Matrix;
    shift_type::Vector{String},
    density_damp::Float64,
    wave_damp::Float64,
    n_genes::Int64)
    
    gene_set = collect(1:n_genes)
    for type in shift_type
        local_type = getfield(base_cells, Symbol(type))
        for cell in keys(local_type)
            if type == "type_shifts" && Int(cell_label) ≠ Int(cell)
                continue
            end
            local_noise = getindex(local_type[cell],1).noise
            local_set = getindex(local_type[cell],2)
            local_shift = getindex(local_type[cell],3)
            # if type == "chain_shifts"
                
            # else
                local_cell = ecosystem[cell,:]
            #end
           
            
            for (gene, shift) in zip(local_set, local_shift)
                local_shift = aggregate_shift(shift, local_cell, type, density_damp, wave_damp)
                local_shift = add_noise(local_shift, local_noise)
                gene_set = apply_shift(gene_set, gene, local_shift)
            end
        end
    end
    return gene_set
end



function add_ecosystems(
    tissue::Tissue,
    base_cells::BaseCells)::Tissue
    cell_graph = tissue.cell_graph
    cell_labels = tissue.cell_labels
    wave_diffusion = base_cells.wave_diffusion
    types = base_cells.types
    density_damp = base_cells.density_damp
    wave_damp = base_cells.wave_damp
    n_genes = base_cells.n_genes
    shift_type = base_cells.shift_type
    ecosystems = Vector{Tuple{Matrix,Vector{Float64}}}(undef,length(cell_labels))
    for eco in eachindex(ecosystems)
        local_niche, dist = get_neighborhood_graph(
            cell_graph,
            eco,
            depth = wave_diffusion)

        local_eco = get_ecosystem(
            local_niche,
            dist,
            cell_labels = cell_labels,
            types = types,
            wave_diffusion = wave_diffusion)
        
        local_shift = ecosystem_shifts(
            base_cells,
            cell_labels[eco],
            local_eco,
            shift_type = shift_type,
            density_damp = density_damp,
            wave_damp = wave_damp,
            n_genes = n_genes)
        ecosystems[eco] = (local_eco,local_shift)
    end
    tissue.ecosystems = ecosystems
    return tissue
end

function collapse_ecosystem(
    ecosystem::Matrix{Float64},
    collapse_to::Int64)
    if collapse_to > 0
        return round.(ecosystem ./ collapse_to) .* collapse_to
    else
        return map(x -> x > 0 ? one(x) : zero(x), ecosystem)
    end
    
end

function trim_ecosystem(
    ecosystem::Matrix{Float64},
    depth_range::Int64 = 3)
    # could throw DomainError here I guess?
    if size(ecosystem, 2) < depth_range
        @warn "Depth range exceeds max depth - fall back to max depth"
        depth_range = size(ecosystem, 2)
    end
    if depth_range < 1
        @warn "Depth range cannot be 0 - fall back to direct neighborhood"
    end
    return ecosystem[:,1:depth_range]

end

function add_ecotypes(
    tissue::Tissue;
    collapse_to::Int64 = 0,
    depth_range::Int64 = 3)::Tissue
    ecosystems = tissue.ecosystems  # Vector{Tuple{Matrix, Vector}}
    cell_labels = tissue.cell_labels
    n = length(ecosystems)

    resize!(cell_labels, n)
    fill!(cell_labels, 0)

    # Storage for unique representatives and their hashes
    unique_hashes = UInt64[]
    unique_matrices = AbstractMatrix[]
    matrix_cell_groups = Vector{Vector{Int}}()

    for i in 1:n
        mat, _ = ecosystems[i]
        mat = collapse_ecosystem(mat, collapse_to)
        mat = trim_ecosystem(mat, depth_range)
        h = hash(mat)

        matched = false
        for (j, existing_h) in enumerate(unique_hashes)
            if h == existing_h && mat == unique_matrices[j]
                push!(matrix_cell_groups[j], i)
                matched = true
                break
            end
        end

        if !matched
            push!(unique_hashes, h)
            push!(unique_matrices, mat)
            push!(matrix_cell_groups, [i])
        end
    end

    for (label, cell_indices) in enumerate(matrix_cell_groups)
        for idx in cell_indices
            cell_labels[idx] = label
        end
    end

    tissue.cell_labels = cell_labels
    return tissue
end
