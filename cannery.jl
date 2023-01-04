# https://github.com/sylvaticus/juliatutorial/blob/master/assets/JuMP.ipynb
# Import of the JuMP, GLPK and CSV modules (the latter one just to import the data from a header based table, as in the original trasnport example in GAMS 
#using Pkg
#Pkg.add("CSV")
using JuMP, GLPK, CSV, DataFrames

# Define sets #
#  Sets
#       i   canning plants   / seattle, san-diego /
#       j   markets          / new-york, chicago, topeka / ;
plants  = ["seattle","san_diego"]          # canning plants
markets = ["new_york","chicago","topeka"]  # markets

# Define parameters #
#   Parameters
#       a(i)  capacity of plant i in cases
#         /    seattle     350
#              san-diego   600  /
a = Dict(              # capacity of plant i in cases
  "seattle"   => 350,
  "san_diego" => 600,
)
 
#       b(j)  demand at market j in cases
#         /    new-york    325
#              chicago     300
#              topeka      275  / ;
b = Dict(              # demand at market j in cases
  "new_york"  => 325,
  "chicago"   => 300,
  "topeka"    => 275,
)

# Table d(i,j)  distance in thousands of miles
#                    new-york       chicago      topeka
#      seattle          2.5           1.7          1.8
#      san-diego        2.5           1.8          1.4  ;
d_table = CSV.read(IOBuffer("""
plants    new_york chicago topeka
seattle   2.5      1.7     1.8
san_diego 2.5      1.8     1.4
"""), DataFrame; delim=" ", ignorerepeated=true, copycols=true)

d = Dict( (r[:plants],m) => r[Symbol(m)] for r in eachrow(d_table), m in markets)
# Here we are converting the table in a "(plant, market) => distance" dictionary
# r[:plants]:   the first key, row field using a fixed header
# m:            the second key
# r[Symbol(m)]: the value, the row field with a dynamic header
 
# Scalar f  freight in dollars per case per thousand miles  /90/ ;
f = 90 # freight in dollars per case per thousand miles 
 
# Parameter c(i,j)  transport cost in thousands of dollars per case ;
#            c(i,j) = f * d(i,j) / 1000 ;
# We first declare an empty dictionary and then we fill it with the values
c = Dict() # transport cost in thousands of dollars per case ;
[ c[p,m] = f * d[p,m] / 1000 for p in plants, m in markets]

# Model declaration (transport model)
trmodel = Model(with_optimizer(GLPK.Optimizer,msg_lev=GLPK.MSG_ON))
# or 
# trmodel = Model(with_optimizer(Ipopt.Optimizer,print_level=0)) # would require "using Ipopt" instead of GLPK

## Define variables ##
#  Variables
#       x(i,j)  shipment quantities in cases
#       z       total transportation costs in thousands of dollars ;
#  Positive Variable x ;
@variables trmodel begin
    x[p in plants, m in markets] >= 0 # shipment quantities in cases
end

## Define contrains ##
# supply(i)   observe supply limit at plant i
# supply(i) .. sum (j, x(i,j)) =l= a(i)
# demand(j)   satisfy demand at market j ;  
# demand(j) .. sum(i, x(i,j)) =g= b(j);
@constraints trmodel begin
    supply[p in plants],   # observe supply limit at plant p
        sum(x[p,m] for m in markets)  <=  a[p]
    demand[m in markets],  # satisfy demand at market m
        sum(x[p,m] for p in plants)   >=  b[m]
end

# Objective
@objective trmodel Min begin
    sum(c[p,m]*x[p,m] 
    for p in plants, m in markets)
end

print(trmodel)

optimize!(trmodel)

status = termination_status(trmodel)
if (status == MOI.OPTIMAL) || (status == MOI.LOCALLY_SOLVED)
    println("Objective value: ",objective_value(trmodel))
    println(value.(x))
    println("Shadow prices of supply:")
    [println("$p = $(dual(supply[p]))") for p in plants]
    println("Shadow prices of demand:")
    [println("$m = $(dual(demand[m]))") for m in markets]
else
    println("Model didn't solved")
    println(status)
end

# Expected result:
# obj= 153.675
#['seattle','new-york']   = 50
#['seattle','chicago']    = 300
#['seattle','topeka']     = 0
#['san-diego','new-york'] = 275
#['san-diego','chicago']  = 0
#['san-diego','topeka']   = 275

