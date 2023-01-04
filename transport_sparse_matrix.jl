# https://shuvomoy.github.io/blogs/posts/Solving-transportation-problem-in-Julia-Jump/#a_test_example_transportation_problem


cities =
[
:BANGKOK; :LONDON; :PARIS; :SINGAPORE;  :NEWYORK;  :ISTANBUL;  :DUBAI;  :KUALALUMPUR;  :HONGKONG;  :BARCELONA
]

products   =
[
:smartphone; :tablet; :laptop
]

capacity    = 700

struct Route
  p::Symbol # p stands for product
  o::Symbol # o stands for origin
  d::Symbol # d stands for destination
end

routesExample =
[
Route(:smartphone,:BANGKOK,:SINGAPORE);
Route(:smartphone,:BANGKOK,:NEWYORK);
Route(:smartphone,:BANGKOK,:ISTANBUL);
Route(:smartphone,:BANGKOK,:DUBAI);
]
routesExample[2] # Will give the second route
routesExample[4].d # Will give the demand city of the 4th route

routes =
[
Route(:smartphone,:BANGKOK,:SINGAPORE);
Route(:smartphone,:BANGKOK,:NEWYORK);
Route(:smartphone,:BANGKOK,:ISTANBUL);
Route(:smartphone,:BANGKOK,:DUBAI);
Route(:smartphone,:BANGKOK,:KUALALUMPUR);
Route(:smartphone,:BANGKOK,:HONGKONG);
Route(:smartphone,:BANGKOK,:BARCELONA);
Route(:smartphone,:LONDON,:SINGAPORE);
Route(:smartphone,:LONDON,:NEWYORK);
Route(:smartphone,:LONDON,:ISTANBUL);
Route(:smartphone,:LONDON,:DUBAI);
Route(:smartphone,:LONDON,:KUALALUMPUR);
Route(:smartphone,:LONDON,:HONGKONG);
Route(:smartphone,:LONDON,:BARCELONA);
Route(:smartphone,:PARIS,:SINGAPORE);
Route(:smartphone,:PARIS,:NEWYORK);
Route(:smartphone,:PARIS,:ISTANBUL);
Route(:smartphone,:PARIS,:DUBAI);
Route(:smartphone,:PARIS,:KUALALUMPUR);
Route(:smartphone,:PARIS,:HONGKONG);
Route(:smartphone,:PARIS,:BARCELONA);

Route(:tablet,:BANGKOK,:SINGAPORE);
Route(:tablet,:BANGKOK,:NEWYORK);
Route(:tablet,:BANGKOK,:ISTANBUL);
Route(:tablet,:BANGKOK,:DUBAI);
Route(:tablet,:BANGKOK,:KUALALUMPUR);
Route(:tablet,:BANGKOK,:HONGKONG);
Route(:tablet,:BANGKOK,:BARCELONA);
Route(:tablet,:LONDON,:SINGAPORE);
Route(:tablet,:LONDON,:NEWYORK);
Route(:tablet,:LONDON,:ISTANBUL);
Route(:tablet,:LONDON,:DUBAI);
Route(:tablet,:LONDON,:KUALALUMPUR);
Route(:tablet,:LONDON,:HONGKONG);
Route(:tablet,:LONDON,:BARCELONA);
Route(:tablet,:PARIS,:SINGAPORE);
Route(:tablet,:PARIS,:NEWYORK);
Route(:tablet,:PARIS,:ISTANBUL);
Route(:tablet,:PARIS,:DUBAI);
Route(:tablet,:PARIS,:KUALALUMPUR);
Route(:tablet,:PARIS,:HONGKONG);
Route(:tablet,:PARIS,:BARCELONA);

Route(:laptop,:BANGKOK,:SINGAPORE);
Route(:laptop,:BANGKOK,:NEWYORK);
Route(:laptop,:BANGKOK,:ISTANBUL);
Route(:laptop,:BANGKOK,:DUBAI);
Route(:laptop,:BANGKOK,:KUALALUMPUR);
Route(:laptop,:BANGKOK,:HONGKONG);
Route(:laptop,:BANGKOK,:BARCELONA);
Route(:laptop,:LONDON,:SINGAPORE);
Route(:laptop,:LONDON,:NEWYORK);
Route(:laptop,:LONDON,:ISTANBUL);
Route(:laptop,:LONDON,:DUBAI);
Route(:laptop,:LONDON,:KUALALUMPUR);
Route(:laptop,:LONDON,:HONGKONG);
Route(:laptop,:LONDON,:BARCELONA);
Route(:laptop,:PARIS,:SINGAPORE);
Route(:laptop,:PARIS,:NEWYORK);
Route(:laptop,:PARIS,:ISTANBUL);
Route(:laptop,:PARIS,:DUBAI);
Route(:laptop,:PARIS,:KUALALUMPUR);
Route(:laptop,:PARIS,:HONGKONG);
Route(:laptop,:PARIS,:BARCELONA);
]

struct Supply
  p::Symbol
  o::Symbol
end

suppliesExample = Supply[] # Creates a 0 element array of immutable type Supply
for r in routesExample # For every element of the route route
  push!(suppliesExample, Supply(r.p, r.o)) # pick the product and origin city and push it in  supplies
end
suppliesExample

costCofExample = [120; 205; 310; 45.0]
costRoutesExample=Dict{Route, Float64}()# Create an empty dictionary
# where the key is Route and the value is Float64
for i in 1:length(routesExample)
  costRoutesExample[routesExample[i]]=costCofExample[i]
  # routesExample[i] is the key, and costCofExample[i] is the value
end
costRoutesExample
costRoutesExample[routesExample[4]]
routesExample[3]
costRoutesExample[Route(:smartphone,:BANGKOK,:ISTANBUL)]

#Creating the array supplies
#                   --------
supplies = Supply[] # Creates a 0 element array of immutable type Supply
for r in routes
  push!(supplies, Supply(r.p, r.o))
end

# Creating suppliedAmount dictionary
#          --------------
#It might be better to create this as a dictionary, where the key is the
# element of the array supplies and the value is the corresponding supplied
#amount


suppliedAmount = Dict{Supply, Float64}()
for s in supplies
  if s.p == :smartphone && s.o == :LONDON
    suppliedAmount[s]=800
  elseif s.p == :smartphone && s.o==:BANGKOK
    suppliedAmount[s]=500
  elseif s.p == :smartphone && s.o==:PARIS
    suppliedAmount[s]=600
  elseif s.p == :tablet && s.o==:BANGKOK
    suppliedAmount[s]=1000
  elseif s.p == :tablet && s.o==:LONDON
    suppliedAmount[s]=1500
  elseif s.p == :tablet && s.o == :PARIS
    suppliedAmount[s]=1700
  elseif s.p == :laptop && s.o == :BANGKOK
    suppliedAmount[s]=150
  elseif s.p == :laptop && s.o == :LONDON
    suppliedAmount[s]=250
  elseif s.p == :laptop && s.o == :PARIS
    suppliedAmount[s]=400
  end #if
end #for

struct Customer
  p::Symbol
  d::Symbol
end

# Creating customers array, which is an array of custom immutable Customer
#          ---------
customers = Customer[]
for r in routes
  push!(customers, Customer(r.p, r.d))
end


demandedAmount = Dict{Customer, Float64}()
for c in customers
  #1
  if c.p==:smartphone && c.d==:SINGAPORE
    demandedAmount[c]=400
    #2
  elseif c.p==:tablet && c.d==:SINGAPORE
    demandedAmount[c]=600
    #3
  elseif c.p==:laptop && c.d==:SINGAPORE
    demandedAmount[c]=90
    #4
  elseif c.p==:smartphone && c.d==:NEWYORK
    demandedAmount[c]=200
    #5
  elseif c.p==:tablet && c.d==:NEWYORK
    demandedAmount[c]=650
    #6
  elseif c.p==:laptop && c.d==:NEWYORK
    demandedAmount[c]=110
    #7
  elseif c.p==:smartphone && c.d==:ISTANBUL
    demandedAmount[c]=100
    #8
  elseif c.p==:tablet && c.d==:ISTANBUL
    demandedAmount[c]=300
    #9
  elseif c.p==:laptop && c.d==:ISTANBUL
    demandedAmount[c]=0
    #10
  elseif c.p==:smartphone && c.d==:DUBAI
    demandedAmount[c]=175
    #11
  elseif c.p==:tablet && c.d==:DUBAI
    demandedAmount[c]=350
    #12
  elseif c.p==:laptop && c.d==:DUBAI
    demandedAmount[c]=65
    #13
  elseif c.p==:smartphone && c.d==:KUALALUMPUR
    demandedAmount[c]=550
    #14
  elseif c.p==:tablet && c.d==:KUALALUMPUR
    demandedAmount[c]=950
    #15
  elseif c.p==:laptop && c.d==:KUALALUMPUR
    demandedAmount[c]=185
    #16
  elseif c.p==:smartphone && c.d==:HONGKONG
    demandedAmount[c]=200
    #17
  elseif c.p==:tablet && c.d==:HONGKONG
    demandedAmount[c]=750
    #18
  elseif c.p==:laptop && c.d==:HONGKONG
    demandedAmount[c]=150
    #19
  elseif c.p==:smartphone && c.d==:BARCELONA
    demandedAmount[c]=275
    #20
  elseif c.p==:tablet && c.d==:BARCELONA
    demandedAmount[c]=600
    #21
  elseif c.p==:laptop && c.d==:BARCELONA
    demandedAmount[c]=200
  end
end


costCof =
[34; 7; 8; 10; 11; 74; 9; 18; 5; 15; 6; 23; 81; 18; 20; 10; 9;
13; 25; 85; 13; 40; 17; 7; 16; 20; 80; 9; 24; 5; 15; 11; 23;
90; 22; 19; 15; 16; 15; 24; 100; 21; 37; 12; 9; 16; 14;
88; 9; 28; 13; 17; 8; 32; 100; 18; 28; 15; 18; 16; 30; 102; 15]

# Creating costRoutes dictionary which contains the costs of the relevant routes
costRoutes=Dict{Route, Float64}()
for i in 1:length(routes)
  costRoutes[routes[i]]=costCof[i]
end

# Creating orig, which takes the product as the input and gives the set of origins of that product

orig = Dict{Symbol, Array}()
for i in 1:length(products)
  dummy_array = Symbol[]
  for j in 1:length(routes)
    #println(i, j, products[i] == routes[j].p)
    if products[i] == routes[j].p
      push!(dummy_array, routes[j].o)
      #println(orig[products[i]])
    else
      #println("Oops, something is not right")
    end #if
  end #for
  orig[products[i]]=unique(dummy_array)
end #for

# Creating dest, which takes the product as the input and gives the set of destinations of that product

dest = Dict{Symbol, Array}()
for i in 1:length(products)
  dummy_array = Symbol[]
  for j in 1:length(routes)
    #println(i, j, products[i] == routes[j].p)
    if products[i] == routes[j].p
      push!(dummy_array, routes[j].d)
      #println(orig[products[i]])
    else
      #println("Oops, something is not right")
    end #if
  end #for
  dest[products[i]]=unique(dummy_array)
end #for

# Load packages
using JuMP, GLPK

# Model name
transpModel = Model(GLPK.Optimizer)

# Variable
@variable(transpModel, opt_prod[routes] >= 0)

# Objective
@objective(transpModel, Min, sum(costRoutes[l]*opt_prod[l] for l in routes))

# First Constraint
for pr in products
  for org in orig[pr]
    @constraint(transpModel, sum(opt_prod[Route(pr, org, de)] for de in dest[pr])
    == suppliedAmount[Supply(pr,org)])
  end
end

#Second Constraint
for pr in products
  for de in dest[pr]
    @constraint(transpModel, sum(opt_prod[Route(pr, org, de)] for org in orig[pr])
    == demandedAmount[Customer(pr,de)])
  end
end

# Final constraint:
for org in cities
  for de in cities
    if org!=de
      @constraint(transpModel,
      sum(
      opt_prod[r] for r in routes
      if r.o == org && r.d==de # This will be used as an filtering condition
        )
        <= capacity)
      else
        continue
      end
    end
  end

  statusMipModel = optimize!(transpModel) # solves the model

  println("The optimal solution is, trans= \n", objective_value(transpModel))

  
