include("distfcm.jl")

function find_nodes_ports_or_le(port_objects, lumped_elements_objects, nodi_coord, escal)
    input_positions = []
    output_positions = []

    for port_object in port_objects
        #@assert length(port_object["inputElement"]) == 3
        ipos = zeros(3)
        ipos[1] = port_object["inputElement"][1] * escal
        ipos[2] = port_object["inputElement"][2] * escal
        ipos[3] = port_object["inputElement"][3] * escal
        push!(input_positions, ipos)
        #@assert length(port_object["outputElement"]) == 3
        opos = zeros(3)
        opos[1] = port_object["outputElement"][1] * escal
        opos[2] = port_object["outputElement"][2] * escal
        opos[3] = port_object["outputElement"][3] * escal
        push!(output_positions, opos)
    end
    ports = Dict(
        :port_start => input_positions,
        :port_end => output_positions,
        :port_nodes => zeros(1,1),
    )
    input_positions_lumped = []
    output_positions_lumped = []
    values = []
    types = []
    R_values = []
    L_values = []
    C_values = []
    N_LUMPED_ELEMENTS = length(lumped_elements_objects)
    if N_LUMPED_ELEMENTS == 0
        lumped_elements = Dict(
            :le_start => [],
            :le_end => [],
            :le_nodes => zeros(1,1),
            :R => [],
            :C => [],
            :L => [],
            :type => []
        )
        #@assert length(input_positions_lumped) == N_LUMPED_ELEMENTS && length(output_positions_lumped) == N_LUMPED_ELEMENTS && length(values) == N_LUMPED_ELEMENTS && length(types) == N_LUMPED_ELEMENTS
    else
        for lumped_element_object in lumped_elements_objects
            #@assert length(lumped_element_object["inputElement"]) == 3
            ipos = zeros(3)
            ipos[1] = lumped_element_object["inputElement"][1] * escal
            ipos[2] = lumped_element_object["inputElement"][2] * escal
            ipos[3] = lumped_element_object["inputElement"][3] * escal
            push!(input_positions_lumped, ipos)
            #@assert length(lumped_element_object["outputElement"]) == 3
            opos = zeros(3)
            opos[1] = lumped_element_object["outputElement"][1] * escal
            opos[2] = lumped_element_object["outputElement"][2] * escal
            opos[3] = lumped_element_object["outputElement"][3] * escal
            push!(output_positions_lumped, opos)

            push!(types, lumped_element_object["type"])

            push!(R_values, haskey(lumped_element_object["rlcParams"], "resistance") ? lumped_element_object["rlcParams"]["resistance"] : 0.0)

            push!(L_values, haskey(lumped_element_object["rlcParams"], "inductance") ? lumped_element_object["rlcParams"]["inductance"] : 0.0)

            push!(C_values, haskey(lumped_element_object["rlcParams"], "capacitance") ? lumped_element_object["rlcParams"]["capacitance"] : 0.0)
        end

        #@assert length(input_positions_lumped) == N_LUMPED_ELEMENTS && length(output_positions_lumped) == N_LUMPED_ELEMENTS && length(values) == N_LUMPED_ELEMENTS && length(types) == N_LUMPED_ELEMENTS && length(R_values) == N_LUMPED_ELEMENTS && length(L_values) == N_LUMPED_ELEMENTS && length(C_values) == N_LUMPED_ELEMENTS

        lumped_elements = Dict(
            :le_start => input_positions_lumped,
            :le_end => output_positions_lumped,
            :le_nodes => zeros(1,1),
            :R => R_values,
            :C => C_values,
            :L => L_values,
            :type => types
        )
    end
    Np = size(ports[:port_start], 1)

    # Initialize the port nodes
    ports[:port_nodes] = zeros(Np, 2)

    # Find port start and end nodes by scaling and calling nodes_find_rev
    ports[:port_nodes][:, 1] = nodes_find_rev(ports[:port_start], nodi_coord)
    ports[:port_nodes][:, 2] = nodes_find_rev(ports[:port_end], nodi_coord)

    Nle = size(lumped_elements[:le_start], 1)

    # Initialize the lumped element nodes
    lumped_elements[:le_nodes] = zeros(Nle, 2)

    # Find lumped element start and end nodes by scaling and calling nodes_find_rev
    lumped_elements[:le_nodes][:, 1] = nodes_find_rev(lumped_elements[:le_start], nodi_coord)
    lumped_elements[:le_nodes][:, 2] = nodes_find_rev(lumped_elements[:le_end], nodi_coord)

    return ports, lumped_elements
end

function nodes_find_rev(Nodes_inp_coord, nodi_centri)
    nodes = zeros(size(Nodes_inp_coord, 1), 1)
    
    for k in range(1,size(Nodes_inp_coord, 1))
        # Compute the distances using distfcm, find minimum distances and corresponding node index
        dist = distfcm(Nodes_inp_coord[k, :], nodi_centri)
        nodes_app = findfirst(vec(dist) .== minimum(dist))
        nodes[k] = nodes_app
    end
    
    return nodes
end
