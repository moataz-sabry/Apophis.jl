using AlgebraOfGraphics, CairoMakie

const publication_themes = Theme(
    fontsize=16, font="CMU Serif",
    Axis=(xlabelsize=20, xgridstyle=:dash, ygridstyle=:dash,
        xtickalign=1, ytickalign=1, yticksize=10, xticksize=10,
        xlabelpadding=-5, xlabel="x", ylabel="y"),
    Legend=(framecolor=(:black, 0.5), bgcolor=(:white, 0.5)),
    Colorbar=(ticksize=16, tickalign=1, spinewidth=0.5),
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