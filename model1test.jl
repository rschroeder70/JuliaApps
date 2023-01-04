#using Pkg
#Pkg.add("DataFrames")
import DataFrames

#Pkg.add("XLSX")
import XLSX

const DATA_DIR = joinpath(@__DIR__, "Data")

excel_df = DataFrames.DataFrame(
    XLSX.readtable(joinpath(DATA_DIR, "ship_to_dist.xlsx"), "Sheet 1")...
)

#= Some attributes of the data frame
DataFrames.size(excel_df)
DataFrames.describe(excel_df)
DataFrames.names(excel_df)
eltype.(excel_df)
=#

# remove EX rows and PR
filter!(row -> row."de.State" != "EX", excel_df)
filter!(row -> row."de.State" != "PR", excel_df)
filter!(row -> row."de.State" != "HI", excel_df)
filter!(row -> !((row."de.City" == "GNANGARA") & (row."de.State" == "WA")), excel_df)

# remove or.State
#DataFrames.rename!(excel_df, "or.State" => "or_state", "de.City" => "de_city", "de.State" => "de_state")
#excel_df[:, Not("or.State")] #Cant get this to work!
DataFrames.select!(excel_df, "or.City", "de.City", "de.State", "Miles")

excel_df.Dest = string.(excel_df."de.City", ", ", excel_df."de.State")
DataFrames.rename!(excel_df, "or.City" => "Orig")

DataFrames.select!(excel_df, :Orig, :Dest, :Miles)

# Change to wide format
excel_df = DataFrames.unstack(excel_df, :Orig, :Miles)
#XLSX.writetable("temp.xlsx", collect(DataFrames.eachcol(excel_df)), DataFrames.names(excel_df))

#Create Vector of the Destinations
DEST = excel_df.Dest
#eltype.(Dest)

ORIG = DataFrames.names(excel_df)
ORIG = ORIG[2:6]

demand = repeat([100.0], 864)

dist = Matrix(excel_df[:, 2:6])
parse.(Float64, string.(dist))

dist = dist' # transpose

using JuMP
import Test
import GLPK

model = Model(GLPK.Optimizer)
	@variable(model, trans[1:length(ORIG), 1:length(DEST)] >= 0)
	@objective(
		model,
		Min,
		sum(
			dist[i, j] * trans[i, j]
			for i in 1:length(ORIG), j in 1:length(DEST)
		)
	)
@constraints(model, begin
	#[i in 1:length(ORIG)], sum(trans[i, :]) == supply[i]
	[j in 1:length(DEST)], sum(trans[:, j]) == demand[j]
end)

optimize!(model)
Test.@test termination_status(model) == MOI.OPTIMAL
Test.@test primal_status(model) == MOI.FEASIBLE_POINT
println("the objective is: ")
println(objective_value(model))
println("The optimal solution is: ")
println(value.(trans))