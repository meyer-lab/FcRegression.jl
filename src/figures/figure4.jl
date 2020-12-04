function figure4()
    ITPFit, odf, effects, ActI_df = regressionResult("ITP"; L0 = 1e-9, f = 4, murine = true)
    melFit, odf, effects, ActI_df = regressionResult("melanoma"; L0 = 1e-9, f = 4, murine = true)
    
    setGadflyTheme()
    fig_Mel_Dep = figureW("melanoma"; IgGX = 2, IgGY = 4, L0 = 1e-9, murine = true, legend = true)
    fig_ITP_Dep = figureW("ITP"; IgGX = 2, IgGY = 4, L0 = 1e-9, murine = true, legend = true)
    Act = plotDepletionSynergy(2, 4; dataType = "ITP", fit = ITPFit, L0 = 1e-9, murine = true, Cellidx = 2)
    Rb = plotDepletionSynergy(2, 4; dataType = "ITP", fit = ITPFit, L0 = 1e-9, murine = true, Cellidx = 2, Rbound = true)
    Fc2 = plotDepletionSynergy(2, 4; dataType = "ITP", fit = ITPFit, L0 = 1e-9, murine = true, Cellidx = 2, Recepidx = 3)
    Fc2_Rb = plotDepletionSynergy(2, 4; dataType = "ITP", fit = ITPFit, L0 = 1e-9, murine = true, Cellidx = 2, Recepidx = 3, Rbound = true)

    draw(SVG("figure4.svg", 14inch, 8inch), plotGrid((2, 4), [fig_ITP_Dep[4] fig_ITP_Dep[6] fig_Mel_Dep[4] fig_Mel_Dep[6] Act Rb Fc2 Fc2_Rb]; widths = [3 4 3 4; 1 1 1 1]))
end
