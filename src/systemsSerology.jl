using DataFrames
using CSV
import StatsPlots


""" Plot HIV1.p66 data in terms of FcgRIIIa vs. ADCC data"""
function plotHIV1p66()
    dataDir = joinpath(dirname(pathof(FcgR)), "..", "data")
    dfF = CSV.read(joinpath(dataDir, "alter-MSB", "data-function.csv"))

    df = FcgR.importAlterMSG()
    dfL = df[occursin.("HIV1.p66", df.Sig), :]
    
    # Gather all of the ADCC data from "data-function.csv"
    newdfF = dfF[dfF[!, :ADCC].!= "NA", :]
    ADCClist = [parse(Float64,x) for x in newdfF[!, :ADCC]] # All ADCC data in one array of ints
    dfADCC = newdfF[:, [:Column1, :ADCC]]
    rename!(dfADCC, [:Subject, :ADCC])
    
    # Note what subjects correspond to the ADCC data
    Subjects = dfADCC[!, :Subject]
    
    #Find corresponding data for FcgRIIIa...HIV1.p66:
  
    #find all rows that contain data for FcgRIIIa.F158.6H.HIV1.p66 & FcgR.IIIa.V158.6H.HIV1.p66, 
    #and create new subsetted dataframes    
    newdfL = dfL[occursin.("FcgRIIIa", dfL.Rec), :]
    newdfLF = newdfL[occursin.("F158", newdfL.Vir), :]
    newdfLV = newdfL[occursin.("V158", newdfL.Vir), :]
    
    #do not use subjects in FcgrIIIa table who had no corresponding ADCC data
    F = newdfLF[(newdfLF.Subject .!= 615167) .& (newdfLF.Subject .!= 930173), :]
    V = newdfLV[(newdfLV.Subject .!= 615167) .& (newdfLV.Subject .!= 930173), :]
    
    #Create combined dataframe
    final = DataFrame(Subject = Subjects, ADCC = ADCClist, F158 = F[!, :Value], V158 = V[!, :Value])
    
    #using StatsPlots
    @StatsPlots.df StatsPlots.final StatsPlots.scatter([:F158 :V158], :ADCC, xlabel = "FcgRIIIA Value", ylabel = "ADCC", title = "HIV1.p66")
end
