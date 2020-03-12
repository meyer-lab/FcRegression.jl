using Pkg
Pkg.add("DataFrames")
using DataFrames
Pkg.add("CSV")
using CSV
Pkg.instantiate()
using FcgR

""" Subset systems serology dataset for HIV1.p66 """
function HIV1p66sub()
    df = FcgR.importAlterMSG()
    
    A = Int64[]
    #find all rows that contain data for HIV1.p66
    for i = 1:size(df, 1)
        if (occursin("HIV1.p66", df.Sig[i]) == true)
            push!(A, i)
        end
    end
    
    #create a subsetted dataframe for HIV1.p66 data
    newdf = df[A, :]
    return newdf
    
    #Note: tried using "newdf = df[in.(df.Sig, occursin("HIV1.p66", df[!, :Sig]), :)]" but has 
        #issues with accessing the strings using df[!, :Sig] and other column syntaxes. This might be an easier way to
        #make this function, but I'm not sure how to access the strings in each row for column Sig without manually iterating through.
end   

Pkg.add("Plots")
Pkg.add("StatsPlots")
using StatsPlots


""" Plot HIV1.p66 data in terms of FcgRIIIa vs. ADCC data"""
function plotHIV1p66()
    dataDir = joinpath(dirname(pathof(FcgR)), "..", "data")
    dfF = CSV.read(joinpath(dataDir, "alter-MSB", "data-function.csv"))
    dfL = HIV1p66sub()
    
    #Gather all of the ADCC data from "data-function.csv"
    newdfF = dfF[dfF[!, :ADCC].!= "NA", :]
    ADCClist = [parse(Float64,x) for x in newdfF[!, :ADCC]] #all ADCC data in one array of ints
    dfADCC = newdfF[:, [:Column1, :ADCC]]
    rename!(dfADCC, [:Subject, :ADCC])
    
    #note what subjects correspond to the ADCC data
    Subjects = dfADCC[!, :Subject]
    
    #Find corresponding data for FcgRIIIa...HIV1.p66:
  
    #find all rows that contain data for FcgRIIIa.F158.6H.HIV1.p66 & FcgR.IIIa.V158.6H.HIV1.p66, 
    #and create new subsetted dataframes
    A = Int64[]
    B = Int64[]
    
    for i = 1:size(dfL, 1)
        if (occursin("FcgRIIIa", dfL.Rec[i]) == true)
            if (occursin("F158", dfL.Vir[i]) == true)
                push!(A, i)
            elseif (occursin("V158", dfL.Vir[i]) == true)
                push!(B, i)
            end
        end
    end
    newdfLF = dfL[A, :]
    newdfLV = dfL[B, :]
    
    #do not use subjects in FcgrIIIa table who had no corresponding ADCC data
    F = newdfLF[(newdfLF.Subject .!= 615167) .& (newdfLF.Subject .!= 930173), :]
    V = newdfLV[(newdfLV.Subject .!= 615167) .& (newdfLV.Subject .!= 930173), :]
    
    #check that all of the subject data lines up
    # can be deleted, but good for a sanity check
    a = 0
    b = 0
    for i = 1:size(Subjects, 1)
        if (Subjects[i] == V.Subject[i])
            a += 1
        end
        if (Subjects[i] == F.Subject[i])
            b += 1
        end 
    end
    @assert a == 179
    @assert b == 179
    
    #Create combined dataframe
    final = DataFrame(Subject = Subjects, ADCC = ADCClist, F158 = F[!, :Value], V158 = V[!, :Value])
    
    #using StatsPlots
    @df final scatter([:F158 :V158], :ADCC, xlabel = "FcgRIIIA Value", ylabel = "ADCC", title = "HIV1.p66")
end

