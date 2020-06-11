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
    df1 = FcgR.antigenTables("6H.HIV1.p66")
    final = rename!(df1, :FcgRIIIa_F158 => :F158, :FcgRIIIa_V158 => :V158)

    #using StatsPlots - New Note: Assuming Gadfly plots can fix this? The ADCC values will plot in list order(not in increasing 
    #order) along y axis with this simplification of the function
    # Reimplement in Gadfly 
    # @df final scatter([:F158 :V158], :ADCC, xlabel = "FcgRIIIA Value", ylabel = "ADCC", title = "HIV1.p66")
end


""" For a given antigen, create a dataframe with Receptor values and ADCC values for each patient """
function antigenTables(s::String)
    dfMSG = FcgR.importLuminex()
    rename!(dfMSG, Dict(:ColNames => "Fc"))

    # Find all rows that contain data for given antigen and create a subsetted dataframe
    df = dfMSG[occursin.(s, dfMSG.Fc), :]

    # Need to differentiate between antigens such "gp120.BAL" and "gp120.BAL.Kif" which will both be found with occursin
    for i = size(df, 1):-1:1
        m = match(Regex(s), df.Fc[i], length(df.Fc[i]) - length(s) + 1)  #match only after index where the antigen name must start
        if m === nothing
            deleterows!(df, i)  #delete all the rows that are not for this specific antigen
        end
    end

    df.Fc = replace.(df.Fc, s => "")  #remove antigen name so we only have rec name for each col
    df.Fc = strip.(df.Fc, ['.'])  #remove trailing "."
    df.Fc = replace.(df.Fc, "." => "_")  #Receptor types now use underscores (e.g. "FcgRIIa_R131" and "FcgRIIa_H131")

    rec = unstack(df, :Subject, :Fc, :Value) #stack all the receptors as columns

    #gather ADCC data
    dataDir = joinpath(dirname(pathof(FcgR)), "..", "data")
    dfF = CSV.read(joinpath(dataDir, "alter-MSB", "data-function.csv"))
    ADCConly = dfF[:, [:Column1, :ADCC]] #only want ADCC data
    ADCC = rename!(ADCConly, [:Subject, :ADCC]) #Subjects were called "Column1" before

    #join ADCC data with antigen receptor data
    allAntigen = join(rec, ADCC, on = :Subject, kind = :inner)

    return allAntigen
end

""" Assemble the Cube """
function createCube()
    dataDir = joinpath(dirname(pathof(FcgR)), "..", "data")
    dfMA = CSV.read(joinpath(dataDir, "alter-MSB", "meta-antigens.csv"))
    dfMS = CSV.read(joinpath(dataDir, "alter-MSB", "meta-subjects.csv"))
    dfMD = CSV.read(joinpath(dataDir, "alter-MSB", "meta-detections.csv"))
    dfMD.detection = replace.(dfMD.detection, "." => "_")
    dfMD.Number = 1:22 #Add a column of index numbers that correspond to each receptor/detection

    Cube = Array{Union{Nothing, Float64}}(nothing, 181, 22, 41) #create a the Cube, filled with nothing
    #Subjects down (181)
    #detections/receptors across (22)
    #antigens each slice (41)

    Subjects = dfMS.Column1

    #Massive for loop that will find correct index for each data point in antigen tables and put into correct index in the Cube

    for p = 1:size(dfMA, 1)
        A = FcgR.antigenTables(dfMA.antigen[p])    #focus on one antigen at a time (one slice of cube)
        B = describe(A)                            #want column names in a listed table for later
        for j = 1:size(dfMD, 1)                    #run through all possible detections/receptors
            for i = 1:size(B, 1)
                if (Symbol(dfMD.detection[j]) == B.variable[i])    #see if current receptor binds to current antigen (exists in table)
                    df2 = select(A, Symbol(dfMD.detection[j]), :Subject)  #create a subsetted dataframe for this receptor antigen combo
                    rename!(df2, [:detection, :Subject])           #need uniform names 
                    for l = 1:size(df2, 1)           #will now match subjects in cube to subjects in dataframe
                        for n = 1:size(Subjects, 1)
                            if (df2.Subject[l] == Subjects[n]) #find index where subject, receptor, and antigen data line up
                                Cube[n, j, p] = df2.detection[l]  #input data point
                            end
                        end
                    end
                end
            end
        end
    end

    return Cube
end
