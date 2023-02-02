using JuMP
using GLPK
using PyPlot



function Run()

n = 3
m = 3

times = [[1, 2, 2],
         [2, 1, 2],
         [1, 2, 1]]
         
         
         
         machines = [[3, 1, 2],
            [2, 3, 1],
            [3, 2, 1]]
            
model=Model(GLPK.Optimizer)      
@variable(model,c >= 0)       
@variable(model, y[1:n, 1:n, 1:m], Bin);
@variable(model, x[1:n, 1:m] >= 0);                
                





M = sum(times[i][j] for i=1:n for j=1:m)

for (j, i) in Iterators.product(1:n, 2:m)
 @constraint(model,x[j,machines[j][i]] - x[j,machines[j][i-1]] >= times[j][machines[j][i-1]])
 end
 
 
 for (j, k) in Iterators.product(1:n, 1:n)
 if k != j
 for i in 1:m
 @constraint(model,x[j,i] - x[k,i] + M*y[j,k,i] >= times[k][i])
 @constraint(model,-x[j,i] + x[k,i] - M*y[j,k,i] >= times[j][i] - M)
 end
 end
 end
 
 
 for j in 1:n
     @constraint(model, c - x[j,machines[j][m]] >= times[j][machines[j][m ]])
     end
     
      @objective(model, Min, c)
 print(model)
     optimize!(model)
     
     
                   
  for (j, i) in Iterators.product(1:n, 1:m)
     println(string("task ",j," starts on machine ",machines[j][i] ," at time", value(x[j,machines[j][i]]), " with time ", times[j][machines[j][i]]))
    end

    println("--------")

  for (j, i) in Iterators.product(1:n, 1:m)
         println(string("task ",j," starts on machine ",i ," at time", value(x[j,i]), " with time ", times[j][i]))
     end

figure()

for (j, i) in Iterators.product(1:n, 1:m)
    a=length(value(x[j,i]).+0.1:0.1:value(x[j,i])+times[j][i].-0.1)

    plot(value(x[j,i]).+0.1:0.1:value(x[j,i])+times[j][i].-0.1,i.*ones(a))
   end



                          
end

Run()

