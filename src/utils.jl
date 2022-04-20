#= to-do list
1. not to create auxillary matrices if the keyword is not to be found at all in the data file
2. better regexes
3. cut down operations for upcoming calculations, ex. W⁻¹
4. offset arrays
5. indexin function?
=#

searchdir(path, key) = filter(x -> contains(x, key), readdir(path)) ### for easily find data files

function isreaction(line) ## checks if the line contains a reaction, used later to assign auxillary parameters to the corresponding reactions
    occursin(r"<?=>?", line)
end

function readfile(file_path) ## filters the data file from comment lines

    file_data = String[]
    for line in eachline(file_path, keep=true)
        if occursin(r"^((?!\s*!))", line)
            commentsout = replace(line, r"!.*" => "")
            push!(file_data, commentsout)
        end
    end
    return file_data
end

function parsemoles(molespecies) ## returns the number of moles, parsed as "Float64", and the species of the reactants and products of a given reaction

    T = Float64
    if occursin(r"^\d", molespecies)
        splitcompound = split(molespecies, r"(?<=^\d)") ## not working with double digits or decimal moles
        molestring = first(splitcompound)
        moles = parse(T, molestring)
        species = last(splitcompound)
    else
        moles = one(T)
        species = molespecies
    end
    return moles, species
end

function findindex(item, itemlist) ## assigns the given species/reaction to an index in the species/reactions vector

    itemindex = findfirst(itemlist) do x ## indexin?
        occursin(item, x)
    end

    return itemindex
end

function findlines(key, list) ## finds the indicies and lines of the given keyword in the given list

    regex = get(regexdictionary, key, nothing) ## indexin?
    keyindices = findall(list) do x
        occursin(regex, x)
    end
    keylines = list[keyindices]
    return keylines, keyindices
end

function decompose(equation) ## decomposes an equation into its reactants and products
    equationsplit = split(equation, r"\s*<?=>?\s*", keepempty=false)
    lhs = first(equationsplit)
    rhs = last(equationsplit)

    regex = r"\(?\+m?\)?|\s+"i
    reactants = split(lhs, regex, keepempty=false)
    products = split(rhs, regex, keepempty=false)
    return reactants, products
end

function findblock(blocktitel, data) ## returns the data of a block from the given data file, example: SPECIES or ELEMENTS block

    regex = Regex("^\\h*$blocktitel", "i")
    ifirst = findindex(regex, data)
    ifirst += 1

    ilast = findnext(data, ifirst) do x
        occursin(r"^end"i, x)
    end
    isnothing(ilast) ? ilast = length(data) - 1 : ilast -= 1

    block = data[ifirst:ilast]
    return block
end

function extractblock(blocktitel, data) ## extracts data out of elements and species blocks

    block = findblock("$blocktitel", data)
    blocklength = length(block)

    if isone(blocklength)
        onlyline = only(block)
        blockvector = split(onlyline)
    else
        joinrange = join(block)
        blockvector = split(joinrange)
    end

    return blockvector
end

function extractequation(reactionline) ## extracts the reaction equation from the reaction line
    splitreactionline =
        split(reactionline, r"(?<!(<|=|>|\+))[ \t]+(?!\g<1>)", keepempty=false)
    equation = first(splitreactionline)
    return equation
end

function extractreactions(type, reactionsblock) ## extracts the reaction lines from the reactions block of type: elementary, threebody or falloff

    reactionslines, _ = findlines(type, reactionsblock)
    reactions::Vector{SubString{String}} = extractequation.(reactionslines) ##otherwise type Any?
    return reactions
end

function indexreaction(auxillaryindex, reactionsblock, assignreactions) ## assigns auxillary parameters to the corresponding reactions

    previndex = findprev(isreaction, reactionsblock, auxillaryindex)
    prevreactionline = reactionsblock[previndex]

    prevreaction = extractequation(prevreactionline)
    assignindex = findindex(prevreaction, assignreactions)
    return assignindex
end

function extractauxillary(keyword, reactionsblock, assignreactions) ## fills a matrix with auxillary parameters

    isequal(keyword, :troe) ? totalparameters = 4 : totalparameters = 3

    auxlines, auxindicies = findlines(keyword, reactionsblock)

    totalreactions = length(assignreactions)
    auxmatrix = zeros(Float64, totalreactions, totalparameters)

    for i in eachindex(auxindicies)
        assignindex = indexreaction(auxindicies[i], reactionsblock, assignreactions)
        auxparameters =
            split(auxlines[i], r"^.*?(?=[\s+-\/]\d(?:\.|\d|\s))|[\s\/]"i, keepempty=false) #r"[ /]|(?i)troe(?-i)"   \s+|[a-df-z\/]|(?<!\d)e

        for parameter in eachindex(auxparameters)
            auxmatrix[assignindex, parameter] = parse(Float64, auxparameters[parameter])
        end
    end
    return auxmatrix
end

function extractalpha(specieslist, reactionsblock, assignreactions) ## fills a matrix with limiting factors

    alphalines, alphaindices = findlines(:alpha, reactionsblock)

    totalspecies = length(specieslist)
    totalreactions = length(assignreactions)
    alphamatrix = ones(Float64, totalreactions, totalspecies)

    for l in eachindex(alphalines)
        reactionindex = indexreaction(alphaindices[l], reactionsblock, assignreactions)
        speciesfactors = split(alphalines[l], r"[\s/]", keepempty=false)
        totalfactors = length(speciesfactors)

        for i = 1:2:totalfactors
            species = speciesfactors[i]
            factorvalue = speciesfactors[i+1]

            speciesindex = findindex(species, specieslist)
            parsevalue = parse(Float64, factorvalue)
            alphamatrix[reactionindex, speciesindex] = parsevalue
        end
    end
    return alphamatrix
end

function assignstoichiometry!(stoichiomatrix, reactionindex, compounds, specieslist) ## assigns stoichiometric ratios to the corresponding species in stoichiometry matrix

    for compound in compounds
        moles, species = parsemoles(compound)
        regex = Regex("^\\Q$(species)\\E\$", "i") ## fix this regarding findindex
        #show(regex)
        speciesindex = findindex(regex, specieslist)
        stoichiomatrix[reactionindex, speciesindex] += moles
    end
    return nothing
end

function extractstoichiometry(reactionlist, specieslist) ## fills matrices with stoichiometry

    totalspecies = length(specieslist)
    totalreactions = length(reactionlist)

    stoichioreactants = zeros(Float64, totalreactions, totalspecies)
    stoichioproducts = zeros(Float64, totalreactions, totalspecies)

    for reactionindex in eachindex(reactionlist)
        reactants, products = decompose(reactionlist[reactionindex])
        assignstoichiometry!(stoichioreactants, reactionindex, reactants, specieslist)
        assignstoichiometry!(stoichioproducts, reactionindex, products, specieslist)
    end

    stoichiosmatrix = stoichioproducts - stoichioreactants
    stoichiosum = sum(stoichiosmatrix, dims=2)

    return sparse(stoichiosmatrix), stoichioreactants, stoichioproducts, stoichiosum
end

function stoichiometric_indices(stoichiosmatrix) ## could be integrated in previous function

    A = stoichiosmatrix
    indicies = [Int[] for _ in axes(A, 1)]
    for i in axes(A, 2), k in axes(A, 1)
        iszero(A[k, i]) || push!(indicies[k], i)
    end
    return indicies
end

function extractweight(speciesline) ## computes the molecular weight from the given species line

    moleweight = 0.0
    regex = r"(?:H|O|N|C|AR|HE)\s+\d(?=H|O|N|C|AR|HE|\h)"i ## bad?
    matchpairs = eachmatch(regex, speciesline)

    for pair in matchpairs
        el, moles = split(pair.match)
        elweight = get(elementsweight, uppercase(el), nothing) ## uppercase to avoind Ar not being found
        elmoles = parse(Float64, moles)
        moleweight += elweight * elmoles
    end
    return moleweight
end

function findlineone(species, specieslist, thermoblock) ## finds line number one of a species in the thermodynamics block

    speciesregex = Regex("^\\s*(\\Q$(specieslist[species])\\E)\\s+.*1\\s*\$", "i") ## "^\\h*(\\Q$(specieslist[species])\\E)\\h+.*1\\h*\$" old regex
    speciesindex = findindex(speciesregex, thermoblock)
    speciesline = thermoblock[speciesindex]
    return speciesline, speciesindex
end

function findtemperature(speciesline) ## finds the middle temperature in line 1 for a species in the thermodynamics block
    tempmatch = match(r"\d+\.\d*(?=\s+1\D)", speciesline)
    #show(speciesline)
    temperature = parse(Float64, tempmatch.match)
    return temperature
end

function sortcoeffs(speciesindex, thermoblock) ## sorts the thermodynamic coeffs for better data processing
    joincoefflines = join(thermoblock[speciesindex+1:speciesindex+3])
    splitcoeffs = split(joincoefflines, r"(?<!e)(?=-)|\s+|\d\s*\n"i, keepempty=false)
    #show(splitcoeffs)
    return splitcoeffs
end

function extractthermo(specieslist, thermoblock) ## extracts the thermodynamic coeffs, molecular weight and middle temperature from the given thermo block
    totalspecies = length(specieslist)
    totalcoefficients = 7
    highcoeffmatrix = zeros(Float64, totalcoefficients, totalspecies)
    lowcoeffmatrix = zeros(Float64, totalcoefficients, totalspecies)

    moleweight = zeros(Float64, totalspecies)
    midtemperatures = zeros(Float64, totalspecies)

    for species in eachindex(specieslist)
        speciesline, speciesindex = findlineone(species, specieslist, thermoblock)

        midtemperatures[species] = findtemperature(speciesline)
        moleweight[species] = extractweight(speciesline)

        splitcoeffs = sortcoeffs(speciesindex, thermoblock)

        for coeff in 1:totalcoefficients
            highcoeffmatrix[coeff, species] = parse(Float64, splitcoeffs[coeff])
            lowcoeffmatrix[coeff, species] = parse(Float64, splitcoeffs[coeff+totalcoefficients])
        end
    end
    return moleweight, midtemperatures, highcoeffmatrix, lowcoeffmatrix
end

# function readmec(title; print=false)

#     #title = Symbol(title) ##???
#     #if isdefined(Apophis, title) == false
#     strng = string(title)
#     #expr = :(const $title = mechanism($strng))
#     expr = :(mechanism = mechan($strng))
#     eval(expr)
#     #Base.eval(Main, expr)
#     #end

#     if print == true
#         #exprs = title
#         #exprs = :(Main.$title)
#         #mechanism = eval(exprs)
#         println(title, ":", "\t$(length(mechanism.elements)) elements", "\t$(length(mechanism.species)) species", "\t$(length(mechanism.reactions)) reactions")
#         println("\t----------", "\t----------", "\t-------------")
#         n = 3
#         for i in 1:n
#             println("\t$(mechanism.elements[i])\t", "\t$(mechanism.species[i])\t", "\t$(mechanism.reactions[i])")
#             i == n && println("\t...\t", "\t...\t", "\t...")
#         end
#     end
#     return nothing
# end

function readmechanism(title) ## main function, calls previous functions in order

    #mechanism_path = pkgdir(Apophis, "mechanisms")
    mechanism_path = "C:/Users/m7mdf/Apophis.jl/mechanisms"
    mech_file_path = mechanism_path * "/$(title)_mech.dat"
    thermo_file_path = mechanism_path * "/$(title)_thermo.dat"
    #trans_file_path = mechanism_path * "/$(title)_trans.dat"

    mech_data = readfile(mech_file_path)

    if isfile(thermo_file_path)
        thermo_data = readfile(thermo_file_path)
        thermoblock = findblock("thermo", thermo_data)
    else
        thermoblock = findblock("thermo", mech_data)
    end

    elements = extractblock("elements", mech_data)
    species = extractblock("species", mech_data)

    molecular_weight, common_temperature, upper_temperature_coefficients, lower_temperature_coefficients = extractthermo(species, thermoblock)
    inverse_molecular_weight = inv.(molecular_weight)

    reactionsblock = findblock("reactions", mech_data)

    elementary_reactions = extractreactions(:elementary, reactionsblock)
    threebody_reactions = extractreactions(:threebody, reactionsblock)
    falloff_reactions = extractreactions(:falloff, reactionsblock)

    reactions = vcat(elementary_reactions, threebody_reactions, falloff_reactions)
    isreversible = findall(x -> occursin(r"<=>|(?<!<)=(?!>)", x), reactions)
    reversible_parameters = findall(x -> occursin(r"rev"i, x), reactions) ###optimizable
    reversible_equilibrium = setdiff(isreversible, reversible_parameters)

    stoichiometric_reactions, stoichiometric_reactants, stoichiometric_products, stoichiometric_sum =
        extractstoichiometry(reactions, species)


    stoichiometric_transpose = sparse(stoichiometric_reactions') ## fix
    reactants_indices = findall(x -> !iszero(x), stoichiometric_reactants) ##matrix as sparse? or vector
    products_indices = findall(x -> !iszero(x), stoichiometric_products)

    #reactants_indices = stoichiometric_indices(stoichiometric_reactants)
    #products_indices = stoichiometric_indices(stoichiometric_products)

    high_pressure_parameters = extractauxillary(:high, reactionsblock, reactions)
    reverse_reaction_parameters = extractauxillary(:rev, reactionsblock, reactions)
    low_pressure_parameters = extractauxillary(:low, reactionsblock, falloff_reactions)
    troe_parameters = extractauxillary(:troe, reactionsblock, falloff_reactions)

    alphareactions = vcat(threebody_reactions, falloff_reactions) ## potential reactions with enhancement factors
    enhancement_factors = extractalpha(species, reactionsblock, alphareactions)

    return Mechanism(
        elements,
        species,
        reactions,
        elementary_reactions,
        threebody_reactions,
        falloff_reactions,
        common_temperature,
        molecular_weight,
        inverse_molecular_weight,
        reversible_parameters,
        reversible_equilibrium,
        reactants_indices,
        products_indices,
        lower_temperature_coefficients,
        upper_temperature_coefficients,
        high_pressure_parameters,
        low_pressure_parameters,
        troe_parameters,
        reverse_reaction_parameters,
        enhancement_factors,
        stoichiometric_products,
        stoichiometric_reactants,
        stoichiometric_reactions,
        stoichiometric_transpose,
        stoichiometric_sum
    )
end