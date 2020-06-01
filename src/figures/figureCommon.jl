"""
Whenever making a plot with Gadfly, use `style()` instead of `Theme()` to add effects
"""
function setGadflyTheme()
    Gadfly.push_theme(Theme(
        background_color = nothing,
        grid_line_style = :dot,
        panel_stroke = colorant"black",
        plot_padding = [0.02w, 0.02w, 0.07h, 0.0pt],
        major_label_font = "Helvetica",
        major_label_color = colorant"black",
        minor_label_font = "Helvetica",
        minor_label_color = colorant"black",
        key_title_font = "Helvetica",
        key_title_color = colorant"black",
        key_label_font = "Helvetica",
        key_label_color = colorant"black",
        point_label_font = "Helvetica",
        point_label_color = colorant"black",
    ))
end


function plotGrid(grid_dim = (1, 1), pls = [], ptitle = nothing)
    @assert length(grid_dim) == 2
    nplots = prod(grid_dim)
    if length(pls) != nplots
        @warn "The number of plots doesn't match the dimension of grid"
    end
    grid = Matrix{Union{Plot, Compose.Context}}(fill(context(), grid_dim))
    for (i, pl) in enumerate(pls)
        if i > nplots
            break
        end
        grid[(i - 1) ÷ grid_dim[2] + 1, (i - 1) % grid_dim[2] + 1] = compose(
            context(),
            (context(), text(0.0, 0.0, 'a' - 1 + i, hleft, vtop), font("Helvetica-Bold"), fontsize(30pt), fill(colorant"black")),
            (context(), render(pl)),
        )
    end
    fpl = gridstack(grid)
    if ptitle != nothing
        fpl = title(fpl, ptitle, Compose.font("Helvetica"), Compose.fontsize(20pt), fill(colorant"black"))
    end
    fpl = compose(context(), (context(), fpl), compose(context(), rectangle()), fill(colorant"white"))
    return fpl
end
