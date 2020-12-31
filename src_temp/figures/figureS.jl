function figureS(Cellidx; L0 = 1e-9, f = 4, murine = true)
    setGadflyTheme()
    p1 = plotDepletionSynergy(1, 2; L0 = L0, f = f, murine = murine, Cellidx = Cellidx)
    p2 = plotDepletionSynergy(1, 3; L0 = L0, f = f, murine = murine, Cellidx = Cellidx)
    p3 = plotDepletionSynergy(1, 4; L0 = L0, f = f, murine = murine, Cellidx = Cellidx)
    p4 = plotDepletionSynergy(1, 5; L0 = L0, f = f, murine = murine, Cellidx = Cellidx)
    p5 = plotDepletionSynergy(2, 3; L0 = L0, f = f, murine = murine, Cellidx = Cellidx)
    p6 = plotDepletionSynergy(2, 4; L0 = L0, f = f, murine = murine, Cellidx = Cellidx)
    p7 = plotDepletionSynergy(2, 5; L0 = L0, f = f, murine = murine, Cellidx = Cellidx)
    p8 = plotDepletionSynergy(3, 4; L0 = L0, f = f, murine = murine, Cellidx = Cellidx)
    p9 = plotDepletionSynergy(3, 5; L0 = L0, f = f, murine = murine, Cellidx = Cellidx)
    p10 = plotDepletionSynergy(4, 5; L0 = L0, f = f, murine = murine, Cellidx = Cellidx)

    return p1, p2, p3, p4, p5, p6, p7, p8, p9, p10
end
