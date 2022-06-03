using AlgebraOfGraphics, CairoMakie

# Conservative 7-color palette
# Wong, Bang. "Points of view: Color blindness." (2011): 441.
# https://www.nature.com/articles/nmeth.1618?WT.ec_id=NMETH-201106

paperred = RGB(185 / 255, 22 / 255, 70 / 255)
paperblue = RGB(75 / 255, 101 / 255, 135 / 255)
papergreen = RGB(0 / 255, 158 / 255, 115 / 255)

function papercolors()
    return [
        RGB(185 / 255, 22 / 255, 70 / 255), # red
        RGB(75 / 255, 101 / 255, 135 / 255), # blue
        RGB(0 / 255, 158 / 255, 115 / 255), # green
        RGB(230 / 255, 159 / 255, 0 / 255), # orange
        RGB(204 / 255, 121 / 255, 167 / 255), # reddish purple
        RGB(86 / 255, 180 / 255, 233 / 255), # sky blue
        RGB(213 / 255, 94 / 255, 0 / 255), # vermillion
        RGB(240 / 255, 228 / 255, 66 / 255), # yellow
    ]
end

const publication_themes = Theme(
    fontsize=16, font="CMU Serif",
    Axis=(xlabelsize=20, xgridstyle=:dash, ygridstyle=:dash,
        xtickalign=1, ytickalign=1, yticksize=10, xticksize=10,
        xlabelpadding=-5, xlabel="x", ylabel="y"),
    Legend=(framecolor=(:black, 0.5), bgcolor=(:white, 0.5)),
    Colorbar=(ticksize=16, tickalign=1, spinewidth=0.5),
    palette=(
        color=papercolors(),
        patchcolor=papercolors(),
        marker=[:circle, :utriangle, :cross, :rect, :diamond, :dtriangle, :pentagon, :xcross],
        linestyle=[:solid, :dash, :dot, :dashdot, :dashdotdot],
        side=[:left, :right],
    )
)

CairoMakie.activate!(type="svg")

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

