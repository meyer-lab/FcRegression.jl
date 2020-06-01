function figure2()
    setGadflyTheme()
    fig_melanoma = figureW("melanoma"; IgGX = 2, IgGY = 3, L0 = 1e-9, f = 4, murine = true)
    fig_ITP = figureW("ITP"; IgGX = 1, IgGY = 2, L0 = 1e-8, f = 10, murine = true)
    fig_blood = figureW("blood"; IgGX = 1, IgGY = 2, L0 = 1e-10, f = 4, murine = true)
    fig_bone = figureW("bone"; IgGX = 1, IgGY = 2, L0 = 1e-10, f = 4, murine = true)
    fig_HIV = figureW("HIV"; IgGX = 1, IgGY = 2, L0 = 1e-8, f = 24, murine = true)

    draw(
        SVG("figure2.svg", 20inch, 7inch),
        plotGrid(
            (2, 5),
            [fig_melanoma[1] fig_ITP[1] fig_blood[1] fig_bone[1] fig_HIV[1] fig_melanoma[2] fig_ITP[2] fig_blood[2] fig_bone[2] fig_HIV[2]],
        ),
    )
end

function figure3()
    setGadflyTheme()
    fig_melanoma = figureW("melanoma"; IgGX = 2, IgGY = 3, L0 = 1e-9, f = 4, murine = true)
    fig_ITP = figureW("ITP"; IgGX = 1, IgGY = 2, L0 = 1e-8, f = 10, murine = true)
    fig_blood = figureW("blood"; IgGX = 1, IgGY = 2, L0 = 1e-10, f = 4, murine = true)
    fig_bone = figureW("bone"; IgGX = 1, IgGY = 2, L0 = 1e-10, f = 4, murine = true)
    fig_HIV = figureW("HIV"; IgGX = 1, IgGY = 2, L0 = 1e-8, f = 24, murine = true)

    draw(SVG("figure3.svg", 15inch, 4inch), plotGrid((1, 5), [fig_melanoma[3] fig_ITP[3] fig_blood[3] fig_bone[3] fig_HIV[3]]))
end
