#using AlgebraOfGraphics, CairoMakie

# Conservative 7-color palette
# Wong, Bang. "Points of view: Color blindness." (2011): 441.
# https://www.nature.com/articles/nmeth.1618?WT.ec_id=NMETH-201106

papercolors =
    [RGB(185 / 255, 22 / 255, 70 / 255), # red
        RGB(0 / 255, 158 / 255, 115 / 255), # green
        RGB(75 / 255, 101 / 255, 135 / 255), # blue
        RGB(230 / 255, 159 / 255, 0 / 255), # orange
        RGB(204 / 255, 121 / 255, 167 / 255), # reddish purple
        RGB(86 / 255, 180 / 255, 233 / 255), # sky blue
        RGB(213 / 255, 94 / 255, 0 / 255), # vermillion
        RGB(240 / 255, 228 / 255, 66 / 255), # yellow
        RGBA(247 / 255, 246 / 255, 241 / 255) # white
    ]

const publication_themes = Theme(
    fontsize=16,
    font="Computer Modern",
    colormap=:batlow,
    linecolor=:gray25,
    markercolor=:gray25,
    patchcolor=:gray25,
    Axis=(
        xgridstyle=:dash,
        xgridvisible=false,
        ygridstyle=:dash,
        ygridvisible=false,
        topspinevisible=true,
        topspinecolor=:black, #darkgray
        rightspinevisible=true,
        rightspinecolor=:black, #darkgray
        bottomspinecolor=:black, #darkgray
        leftspinecolor=:black, #darkgray
        xtickalign=1,
        ytickalign=1,
        xticksize=10,
        yticksize=10,
        xtickcolor=:darkgray,
        ytickcolor=:darkgray,
        xticklabelfont="CMU Serif", #light
        yticklabelfont="CMU Serif",
        xlabelfont="CMU Serif", #medium
        ylabelfont="CMU Serif",
        titlefont="CMU Serif"),
    palette=(
        color=papercolors,
        patchcolor=papercolors,
        marker=[:circle, :utriangle, :cross, :rect, :diamond, :dtriangle, :pentagon, :xcross],
        linestyle=[:solid, :dash, :dot, :dashdot, :dashdotdot],
        side=[:left, :right],
    )
);

function forward_profile(solution, index)
    set_aog_theme!()

    t = sol.t
    T = sol[index, :]
    grp = [fill("a", length(t))]

    df = (; t, T, grp)
    ut = data(df)
    layers = visual(Lines)
    plot = layers * ut * mapping(:t => "time [s]", :T => "temperature [K]")

    return draw(plot)
end

function plotsensitivity(gas, dJdg, nr)
    reactions = @. latexstring(replace(gas.mechanism.reactions[1:nr], r"(?<=[a-z])(?=\d)"i => "_", r"="i => "\\longrightarrow"))
    p = sortperm(dJdg, by=abs, rev=false)

    f = Figure()
    ax = Axis(f[1, 1], xlabel=L"S_A [-]", xticks=0.0:0.25:1.0, yticks=(1:nr, reactions))

    barplot!(ax, last(dJdg[p], nr),
        direction=:x,
        color=papercolors[3], ygridstyle=:solid,
        ygridvisible=true, ygridcolor=:lightgray)
    return f
end

function plotsubspaces(W)
    n = size(W, 1)
    x = collect(1:n)
    f = Figure()
    ax = Axis(f[1, 1], xlabel=L"Reactions", xticks=1:n)
    stem!(ax, x, W[:, end],
        marker=:rect,
        color=:white,
        strokecolor=papercolors[1],
        stemcolor=papercolors[1],
        stemwidth=1.25,
        strokewidth=2,
        trunkcolor=:darkgray,
        trunkwidth=1.0,
        trunklinestyle=:dot)

    stem!(ax, x, W[:, end-1],
        marker=:circle,
        color=:white,
        strokecolor=papercolors[3],
        stemcolor=papercolors[3],
        stemwidth=1.25,
        stemlinestyle=:dash,
        strokewidth=2,
        trunkcolor=:darkgray,
        trunkwidth=1.0,
        trunklinestyle=:dot)
    return f
end

function plotpdf(fmore, G, fLLAMmore)
    f = Figure()
    ax = Axis(f[1, 1], xlabel=L"IDT", ylabel=L"PDF")

    fmoreplot = density!(ax, fmore, color=(papercolors[3], 0.075), linestyle=:solid,
        strokecolor=(papercolors[3], 1.0), strokewidth=1.25, label="Monte Carlo")
    gplot = density!(ax, G, color=(papercolors[1], 0.075), linestyle=:solid,
        strokecolor=(papercolors[1], 1.0), strokewidth=1.25, label="Active Subspaces")
    LLAMplot = density!(ax, fLLAMmore, color=(papercolors[2], 0.075), linestyle=:solid,
        strokecolor=(papercolors[2], 1.0), strokewidth=1.25, label="LLAM")
    f
    return f
end
