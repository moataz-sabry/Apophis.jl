#= to-do list
1. not to create auxillary matrices if the keyword is not to be found at all in the data file
2. better regexes
3. cut down operations for upcoming calculations, ex. W⁻¹
4. offset arrays
5. indexin function?
6. duplicate reactions causes error
7. use filter function
8. something to tell why it could'nt read the mechanism
=#

const regexdictionary = Dict(
    :commentout => r"^((?!\s*!))",
    :elementary => r"(?!.*\+m)(?=.*<?=>?)"i,
    :threebody => r"(?<!\()\+m(?!\))"i,
    :falloff => r"\(\+m\)"i,
    :low => r"low"i,
    :troe => r"troe"i,
    :rev => r"rev"i,
    :alpha => r"^(?!.*(plog|rev|low|troe)).*/.*$"i,
    :high => r"<?=>?"i,
) ## contains regular expressions for extracting data based on a specific keyword

const elementsweight = Dict(
    "H" => 1.00784,
    "O" => 15.9994,
    "C" => 12.0107,
    "N" => 14.0067,
    "AR" => 39.948,
    "HE" => 4.002602,
) ## contains molecular weight for the most common elements


function trapezoid(ƒ::VecOrMat, t::Vector)

    m = size(ƒ, 2)
    I = zeros(m)

    for n in 1:length(t)-1, p = 1:m
        I[p] += 0.5 * (t[n+1] - t[n]) * (ƒ[n+1, p] + ƒ[n, p])
    end

    return m == 1 ? only(I) : I
end

function QoI(sol, tₒ, t∞, tᵣ)

    τ = tₒ - t∞
    ilength = length(last(sol.u))

    t = filter(x -> x > tₒ + 0.01tᵣ, sol.t)
    T = sol(t, idxs=ilength).u

    tₛ = t .- 0.01tᵣ
    Tₛ = sol(tₛ, idxs=ilength).u

    I = (T .- Tₛ) .^ 2

    J = trapezoid(I, t) / 2τ
    return J
end

searchdir(path, key) = filter(x -> contains(x, key), readdir(path)) ### to easily find data files

function filterdata(criteria, data)

    regex = regexdictionary[criteria]
    filtered_data = filter(data) do x
        occursin(regex, x)
    end
    return filtered_data
end
function readfile(file_path) ## filters the data file from comment lines

    file_data = readlines(file_path, keep=true)
    filtered_data = replace.(filterdata(:commentout, file_data), r"!.*" => "") ## quick fix
    return filtered_data
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

    regex = r"\s*\+\s*m(?![a-z])\s*|\s+|\(\+\s*m\)|\+"i #r"\(?\+m?\)?|\s+"i this not working in case +MVOX
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

    block = findblock(blocktitel, data)
    joinblock = ""
    for l in block
        joinblock *= join(l)
    end

    blockvector = split(joinblock)
    return blockvector
end

function extractequation(reactionline) ## extracts the reaction equation from the reaction line
    regex = r"(?<!(<|=|>|\+))[ \t]+(?!\g<1>)" ##maybe split reverse with limit = 4?
    splitreactionline = split(reactionline, regex, keepempty=false)
    equation = first(splitreactionline)
    return equation
end

function extractreactions(type, reactionsblock) ## extracts the reaction lines from the reactions block of type: elementary, threebody or falloff

    reactionslines = filterdata(type, reactionsblock)
    reactions::Vector{SubString{String}} = extractequation.(reactionslines) ##otherwise type Any?
    return reactions
end

function indexreaction(auxillaryindex, reactionsblock, assignreactions) ## assigns auxillary parameters to the corresponding reactions

    previndex = findprev(reactionsblock, auxillaryindex) do x
        occursin(r"<?=>?", x)
    end
    prevreactionline = reactionsblock[previndex]

    prevreaction = extractequation(prevreactionline)
    assignindex = findindex(prevreaction, assignreactions)
    return prevreaction, assignindex
end

function extractauxillary(keyword, reactionsblock, assignreactions) ## fills a matrix with auxillary parameters

    isequal(keyword, :troe) ? totalparameters = 4 : totalparameters = 3

    auxlines, auxindicies = findlines(keyword, reactionsblock)

    totalreactions = length(assignreactions)
    auxmatrix = zeros(Float64, totalreactions, totalparameters)
    isaux = Int64[]
    for i in eachindex(auxindicies)
        prevreaction, assignindex = indexreaction(auxindicies[i], reactionsblock, assignreactions)
        if auxmatrix[assignindex, 1] != 0.0
            assignindex = findnext(assignreactions, assignindex + 1) do x ## united against duplicates
                occursin(prevreaction, x)
            end
        end
        auxparameters =
            split(auxlines[i], r"^.*?(?=[\s+-\/]\d*\.\d*)|[\s\/]"i, keepempty=false) #r"[ /]|(?i)troe(?-i)"   \s+|[a-df-z\/]|(?<!\d)e
        for parameter in eachindex(auxparameters)
            auxmatrix[assignindex, parameter] = try
                parse(Float64, auxparameters[parameter])
            catch
                println(auxparameters)
            end
        end
        push!(isaux, assignindex)
    end
    return isaux, auxmatrix
end

function extractalpha(specieslist, reactionsblock, assignreactions) ## fills a matrix with limiting factors

    alphalines, alphaindices = findlines(:alpha, reactionsblock)

    totalspecies = length(specieslist)
    totalreactions = length(assignreactions)
    alphamatrix = ones(Float64, totalreactions, totalspecies)

    for l in eachindex(alphalines)
        _, reactionindex = indexreaction(alphaindices[l], reactionsblock, assignreactions)
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
        #println(compound)
        moles, species = parsemoles(compound)
        regex = Regex("^\\Q$(species)\\E\$", "i") ## fix this regarding findindex
        #show(regex)
        speciesindex = findindex(regex, specieslist)
        isnothing(speciesindex) && println(reactionindex)
        isnothing(speciesindex) && println(compound)
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
    stoichiosum = dropdims(sum(stoichiosmatrix, dims=2), dims=2)

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
    regex = r"(?:H|O|N|C|AR|HE)\s+\d(?=H|O|N|C|AR|HE|(?=\s*(?:G|L|S|0)))"i ## bad?
    matchpairs = eachmatch(regex, speciesline)

    for pair in matchpairs
        el, moles = split(pair.match)

        elweight = get(elementsweight, uppercase(el), nothing) ## uppercase to avoid Ar not being found
        elmoles = parse(Float64, moles)
        moleweight += elweight * elmoles
    end
    return moleweight
end

function findlineone(species, specieslist, thermoblock) ## finds line number one of a species in the thermodynamics block

    speciesregex = Regex("^\\s*(\\Q$(specieslist[species])\\E)")#\\s+.*1\\s*\$", "i") ## "^\\h*(\\Q$(specieslist[species])\\E)\\h+.*1\\h*\$" old regex
    speciesindex = findindex(speciesregex, thermoblock)
    speciesline = thermoblock[speciesindex]
    return speciesline, speciesindex
end

function findtemperature(speciesline) ## finds the middle temperature in line 1 for a species in the thermodynamics block
    tempmatches = eachmatch(r"\d+\.\d*", speciesline)
    temperatures = collect(tempmatches)
    #lowtemperature = parse(Float64, temperatures[1].match)
    #hightemperature = parse(Float64, temperatures[2].match)
    commontemperature = parse(Float64, temperatures[3].match)
    return commontemperature
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
        #println(speciesline)
        midtemperatures[species] = findtemperature(speciesline)
        moleweight[species] = extractweight(speciesline)

        splitcoeffs = sortcoeffs(speciesindex, thermoblock)
        #println(splitcoeffs)

        for coeff in 1:totalcoefficients
            highcoeffmatrix[coeff, species] = parse(Float64, splitcoeffs[coeff])
            lowcoeffmatrix[coeff, species] = parse(Float64, splitcoeffs[coeff+totalcoefficients])
        end
    end
    return moleweight, midtemperatures, highcoeffmatrix, lowcoeffmatrix
end

function reversibility(keyword, reactionsblock, assignreactions)

    _, auxindicies = findlines(keyword, reactionsblock)

    totalreactions = zero(auxindicies)

    for i in eachindex(auxindicies)
        _, assignindex = indexreaction(auxindicies[i], reactionsblock, assignreactions)
        totalreactions[i] = assignindex
    end
    return totalreactions
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

function readmechanism(title, T) ## main function, calls previous functions in order

    mechanism_path = pkgdir(Apophis, "test/mechanisms")
    mech_file_path = mechanism_path * "/$(title)_mech.dat"
    thermo_file_path = mechanism_path * "/$(title)_thermo.dat"
    #trans_file_path = mechanism_path * "/$(title)_trans.dat"

    mech_data = readfile(mech_file_path)

    if isfile(thermo_file_path)
        thermo_data = readfile(thermo_file_path)
        thermoblock = findblock(:thermo, thermo_data)
    else
        thermoblock = findblock(:thermo, mech_data)
    end

    elements = extractblock(:elements, mech_data)
    species = extractblock(:species, mech_data)

    molecular_weight, common_temperature, upper_temperature_coefficients, lower_temperature_coefficients = extractthermo(species, thermoblock)
    inverse_molecular_weight = inv.(molecular_weight)

    reactionsblock = findblock(:reactions, mech_data)

    elementary_reactions = extractreactions(:elementary, reactionsblock)
    threebody_reactions = extractreactions(:threebody, reactionsblock)
    falloff_reactions = extractreactions(:falloff, reactionsblock)

    reactions = vcat(elementary_reactions, threebody_reactions, falloff_reactions)
    isreversible = findall(x -> occursin(r"<=>|(?<!<)=(?!>)", x), reactions)

    stoichiometric_reactions, stoichiometric_reactants, stoichiometric_products, stoichiometric_sum =
        extractstoichiometry(reactions, species)

    stoichiometric_transpose = sparse(stoichiometric_reactions')

    reactants_indices = stoichiometric_indices(stoichiometric_reactants)
    products_indices = stoichiometric_indices(stoichiometric_products)

    _, forward_reaction_parameters = extractauxillary(:high, reactionsblock, reactions)
    pressure_independant_parameters = forward_reaction_parameters[1:(length(reactions)-length(falloff_reactions)), :] ## optimize
    high_pressure_parameters = forward_reaction_parameters[length(reactions)-length(falloff_reactions)+1:length(reactions), :] ## optimize
    is_reverse_parameters, reverse_reaction_parameters = extractauxillary(:rev, reactionsblock, reactions)
    _, low_pressure_parameters = extractauxillary(:low, reactionsblock, falloff_reactions)
    is_troe_parameters, troe_parameters = extractauxillary(:troe, reactionsblock, falloff_reactions)

    is_reversible_equilibrium = setdiff(isreversible, is_reverse_parameters)

    alphareactions = vcat(threebody_reactions, falloff_reactions) ## potential reactions with enhancement factors
    enhancement_factors = extractalpha(species, reactionsblock, alphareactions)

    return Mechanism{T}(
        elements,
        species,
        reactions,
        elementary_reactions,
        threebody_reactions,
        falloff_reactions,
        common_temperature,
        molecular_weight,
        inverse_molecular_weight,
        is_reverse_parameters,
        is_reversible_equilibrium,
        is_troe_parameters,
        lower_temperature_coefficients,
        upper_temperature_coefficients,
        pressure_independant_parameters,
        high_pressure_parameters,
        low_pressure_parameters,
        troe_parameters,
        reverse_reaction_parameters,
        enhancement_factors,
        stoichiometric_products,
        stoichiometric_reactants,
        stoichiometric_reactions,
        stoichiometric_transpose,
        stoichiometric_sum,
        reactants_indices,
        products_indices
    )
end

function init(mech::Symbol, temperature::K, pressure::K=Pₐ; s...) where {K<:Number}

    strng = string(mech)
    gasexpr = :(gas = Gas{$K}(readmechanism($strng, $K)))
    eval(gasexpr)

    W⁻¹ = gas.mechanism.inverse_molecular_weight

    ## integrate with concentrations? ##
    T = gas.initial.temperature[1] = temperature
    P = gas.initial.pressure = pressure

    Y = gas.initial.mass_fractions
    X, _, _ = gas.initial.molar_fractions
    C, _, _ = gas.initial.molar_concentrations

    if first(s)[1] == :rand
        Y .= first(s)[2]
    else
        species = String.(first.(collect(s)))
        fractions = last.(collect(s))
        ∑fractions = sum(fractions)

        indicies = indexin(species, gas.mechanism.species)
        in(nothing, indicies) && error("one or more of species are not part of $(mech) mechanism")
        ∑fractions ≈ one(K) || @warn "∑Y = $(∑fractions) ≠ 1!"
        for i in eachindex(indicies)
            j = indicies[i]
            Y[j] = fractions[i]
        end
    end

    W̅ = inv(Y ⋅ W⁻¹)
    ρ = gas.current.density = gas.initial.density = P * W̅ / (R * T)

    for i in eachindex(gas.mechanism.species)
        YW⁻¹ = Y[i] * W⁻¹[i]
        C[i] = YW⁻¹ * ρ
        X[i] = YW⁻¹ * W̅
    end

    #gas.current = gas.initial

    return nothing
end
