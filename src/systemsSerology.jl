using DataFrames
using CSV

""" Subset systems serology dataset for HIV1.p66 """
function HIV1p66sub()
    df = FcgR.importAlterMSG()

    # Find all rows that contain data for HIV1.p66
    #create a subsetted dataframe for HIV1.p66 data
    return df[occursin.("HIV1.p66", df.Sig), :]
end   


""" Plot HIV1.p66 data in terms of FcgRIIIa vs. ADCC data"""
function plotHIV1p66()
    df1 = FcgR.antigenTables("HIV1.p66") 
    final = rename!(df1, :FcgRIIIa_F158 => :F158, :FcgRIIIa_V158 =>  :V158)
    
     #using StatsPlots - New Note: Assuming Gadfly plots can fix this? The ADCC values will plot in list order(not in increasing 
                        #order) along y axis with this simplification of the function
     # Reimplement in Gadfly 
     # @df final scatter([:F158 :V158], :ADCC, xlabel = "FcgRIIIA Value", ylabel = "ADCC", title = "HIV1.p66")
end


""" For a given antigen, create a dataframe with Receptor values and ADCC values for each patient """
function antigenTables(string)
    dfMSG = FcgR.importAlterMSG()

    # Find all rows that contain data for given antigen and create a subsetted dataframe
    df = dfMSG[occursin.(string, dfMSG.Sig), :]
    
    #Treat each genetic variant as a different receptor (add the genetic variant to overall receptor name in rec column)
    for i=1:size(df, 1)
        df.Rec[i] = df.Rec[i] * "_" * df.Vir[i]
    end
    
    rec = unstack(df, :Subject, :Rec, :Value) #stack all the receptors as columns
   
    #gather ADCC data
    dataDir = joinpath(dirname(pathof(FcgR)), "..", "data")
    dfF = CSV.read(joinpath(dataDir, "alter-MSB", "data-function.csv"))
    newdfF = dfF[dfF[!, :ADCC].!= "NA", :] #ignore subjects who have no ADCC data
    ADCConly = newdfF[:, [:Column1, :ADCC]] #only want ADCC data
    ADCC = rename!(ADCConly, [:Subject, :ADCC]) #Subjects were called "Column1" before
    
    #join ADCC data with antigen receptor data
    allAntigen = join(rec, ADCC, on = :Subject, kind = :inner)
    
    return allAntigen
end
