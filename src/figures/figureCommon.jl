using Gadfly
using Compose

"""
Whenever making a plot with Gadfly, use `style()` instead of `Theme()` to add effects
"""
function setGadflyTheme()
    Gadfly.push_theme(
        Theme(
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
        ),
    )
end


function plotGrid(grid_dim = (1, 1), pls = [], ptitle = nothing; widths = [], heights = [], sublabels::Union{Bool, Vector{Int}, String} = true)
    @assert length(grid_dim) == 2
    nplots = prod(grid_dim)
    if length(pls) != nplots
        @warn "The number of plots doesn't match the dimension of grid"
    end

    if length(widths) == grid_dim[2]
        widths = repeat(reshape(widths, 1, :), grid_dim[1])
    elseif length(widths) == 0
        widths = ones(grid_dim...)
    end
    widths = widths ./ sum(widths, dims = 2)
    @assert size(widths) == grid_dim "Specified widths, $(size(widths)), should have the same size as the grid, $grid_dim"

    if length(heights) > 0
        @assert length(heights) == grid_dim[1] "Specified heights should have the same size as the grid"
        heights = heights ./ sum(heights)
    else
        heights = ones(grid_dim[1]) ./ grid_dim[1]
    end

    # grid[yi][xi]
    grid = Vector(undef, grid_dim[1])
    for i = 1:grid_dim[1]
        grid[i] = Vector{Union{Plot, Compose.Context}}(fill(context(), grid_dim[2]))
    end

    if sublabels == true
        sublabels = trues(nplots)
    elseif sublabels == false
        sublabels = falses(nplots)
    elseif !(sublabels isa String)
        @assert length(sublabels) == nplots
        sublabels = BitVector(sublabels)
    else
        @warn "In plotGrid(): sublabels cannot be convert to a BitVector"
    end
    letter_label = 'a'

    for i = 1:nplots
        xi = (i - 1) % grid_dim[2] + 1
        yi = (i - 1) ÷ grid_dim[2] + 1
        if i <= length(pls)
            label = ""
            if sublabels isa String
                if 'a' <= sublabels[i] <= 'z'
                    label = sublabels[i]
                end
            elseif sublabels[i]
                label = letter_label
                letter_label += 1
            end
            content = (pls[i] === nothing) ? context() : (context(), render(pls[i]))
            grid[yi][xi] = compose(
                context(0, 0, widths[yi, xi], 1),
                (context(), text(0.0, 0.0, label, hleft, vtop), font("Helvetica-Bold"), fontsize(30pt), fill(colorant"black")),
                content,
            )
        else
            grid[yi][xi] = context(0, 0, widths[yi, xi], 1)
        end
    end

    fplrows = Vector(undef, grid_dim[1])
    for i = 1:grid_dim[1]
        fplrows[i] = compose(context(0, 0, 1, heights[i]), hstack(grid[i]))
    end
    fpl = vstack(fplrows...)

    if ptitle !== nothing
        fpl = title(fpl, ptitle, Compose.font("Helvetica"), Compose.fontsize(20pt), fill(colorant"black"))
    end
    fpl = compose(context(), (context(), fpl), compose(context(), rectangle()), fill(colorant"white"))
    return fpl
end
