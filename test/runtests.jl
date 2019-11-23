using Test
using Profile
using FcgR


L0 = 1.0e-8
KxStar = 1.0e-10
f = 16
Rtot = [100.0, 1000.0, 10.0, 10000.0]
IgGC = [1.0]
Kav = ones((1, 4))*1.0e6
ActI = [1.0, 1.0, -1.0, 0.0]


@testset "Can successfully assemble the parameters." begin
    FcgR.polyfc(L0, KxStar, f, Rtot, IgGC, Kav, ActI)
    @time FcgR.polyfc(L0, KxStar, f, Rtot, IgGC, Kav, ActI)

    for i in 1:50
    	@profile FcgR.polyfc(L0, KxStar, f, Rtot, IgGC, Kav, ActI)
    end

    Profile.print(noisefloor=5.0)
end