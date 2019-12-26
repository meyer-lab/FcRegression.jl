
""" Plot an example isobologram. """
function plotIsobologram()
    Kav = importKav(murine=false)
    FcExpr = zeros(6);
    FcExpr[5] = 1000.0;
    
    output = calculateIsobologram(2, 3, 24, 1.0e-8, FcExpr, Kav)
    
    X = range(0,stop=1,length=length(output))

    plot(X, output, title="Receptor Binding vs IgG Composition", xticks=false, legend=false, dpi=72)
    plot!([0, 1], [output[1], output[33]])
    annotate!([(0, 0, text("100% hIgG2",8,:right, rotation=45)),(1.0, 0, text("100% hIgG3",8,:right, rotation=45))])
    ylabel!("hFcgRIIIA-158V Binding")
    xlabel!("Percent hIgG3")
    ylims!((-1, maximum(output) * 1.1))
end

""" Plot an example isobologram. """
function plotIsobologramTwo()
    Kav = importKav(murine=true)
    FcExpr = [2571.0, 12886.0, 12563.0, 2371.0]
    ActVIn = ones(4)
    ActVIn[2] = -1.0

    output = calculateIsobologram(2, 3, 4, 1.0e-9, FcExpr, Kav, actV=ActVIn)
    output = output / maximum(output)

    X = range(0,stop=1, length=length(output))
    
    plot(X, output, title="Activity vs IgG Composition", xticks=false, legend=false, dpi=72)
    plot!([0, 1], [output[1], output[33]])
    annotate!([(0, 0, text("100% mIgG2a",8,:right, rotation=45)),(1.0, 0, text("100% mIgG2b",8,:right, rotation=45))])
    ylabel!("?cMO? Predicted Activity")
    xlabel!("Percent mIgG2b")
    ylims!((-0.02, maximum(output) * 1.1))
end

"""Figure shows the affect of increasing immune complex concentration on synergies for each IgG combination"""
function PlotSynGraph()
    Kav = importKav(murine=true)
    df = importRtot()
    FcgR = df[:,2] #2 = mean cMO
    IC = exp10.(range(-12, stop=-6, length=20))
    S = zeros((length(IC), 10))

    for (ii, value) in enumerate(IC)
        A = synergyGrid(4, value, FcgR, Kav)
        h = collect(Iterators.flatten(A))
        S[ii, 1:4] = h[1:4]
        S[ii, 5:7] = h[6:8]
        S[ii, 8:9] = h[11:12]
        S[ii, 10] = h[16]
    end

    plot(IC, S, xaxis=:log, title="Effect of Concentration on Synergy", label=["IgG1/2a" "IgG1/2b" "IgG1/3" "IgG2a/2b" "IgG2a/3" "IgG2b/3"], legend=:topleft, dpi=72)
    xlabel!("IC Concentration")
    ylabel!("Synergy")
end

function figureB1()
    l = @layout [a; b c]

    p1 = plot(rand(100, 4))
    p2 = plot(rand(100, 4))
    p3 = plot(rand(100, 4))
    p = plot(p1, p2, p3, layout=l)

    savefig(p, "figureB1.pdf")
end
