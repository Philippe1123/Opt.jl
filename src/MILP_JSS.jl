using JuMP
#using Cbc
using PyPlot
using Gurobi
using XLSX



function Run()
    nameOfFile = "JSP_dataset_FT06.xlsx"

    ENV["GUROBI_HOME"] = "/home/philippe/gurobi912/linux64"
    xf = XLSX.readxlsx(nameOfFile)
    n = xf["Parameters!B2"]
    m = xf["Parameters!B3"]



    times = xf[string("ProcessingTime!B2:", xf[string("Sheet4!A", m + 1)], n + 1)]

    machines = xf[string("MachineSequence!B2:", xf[string("Sheet4!A", m + 1)], n + 1)]



    #n = 10
    #m = 10
    #=

    times = [[3 2 4 3],
            [4 5 2 2],
            [4 4 3 1],
            [3 4 1 1]]

    machines = [[1 4 2 3],
                [2 3 1 4],
                [3 4 2 1],
                [4 1 2 3]]


    times = [[1, 2, 2],
             [2, 1, 2],
             [1, 2, 1]]

             times = [[2, 1, 2],
             [1, 2, 2],
             [1, 2, 1]]


                      machines = [[3, 1, 2],
                [2, 3, 1],
                [3, 2, 1]]






    #FT06
    times = [vec([1 3 6 7 3 6]),
            vec([8 5 10 10 10 4]),
            vec([5 4 8 9 1 7]),
            vec([5 5 5 3 8 9]),
            vec([9 3 5 4 3 1]),
            vec([3 3 9 10 4 1])]


      machines=[vec([3 1 2 4 6 5]),
                vec([2 3 5 6 1 4]),
                vec([3 4 6 1 2 5]),
                vec([2 1 3 4 5 6]),
                vec([3 2 5 6 1 4]),
                vec([2 4 6 1 5 3])] 

    =#



    #FT10
    #=
    times =[  vec([29    78     9    36    49    11    62    56    44    21]),
    vec([43    90    75    11    69    28    46    46    72    30]),
    vec([91    85    39    74    90    10    12    89    45    33]),
    vec([81    95    71    99     9    52    85    98    22    43]),
    vec([14     6    22    61    26    69    21    49    72    53]),
    vec([84     2    52    95    48    72    47    65     6    25]),
    vec([46    37    61    13    32    21    32    89    30    55]),
    vec([31    86    46    74    32    88    19    48    36    79]),
    vec([76    69    76    51    85    11    40    89    26    74]),
    vec([85    13    61     7    64    76    47    52    90    45])]


      machines =  [vec([1  2  3   4  5   6   7   8   9  10]),
      vec([1  3  5  10  4   2   7   6   8   9]),
      vec([2  1  4   3  9   6   8   7  10   5]),
      vec([2  3  1   5  7   9   8   4  10   6]),
      vec([3  1  2   6  4   5   9   8  10   7]),
      vec([3  2  6   4  9  10   1   7   5   8]),
      vec([2  1  4   3  7   6  10   9   8   5]),
      vec([3  1  2   6  5   7   9  10   8   4]),
      vec([1  2  4   6  3  10   7   8   5   9]),
      vec([2  1  3   7  9  10   6   4   5   8]),]

    =#
    #=
    #LA19
     times=[vec([44 5 58 97 9 84 77 96 58 89]),
     vec([15 31 87 57 77 85 81 39 73 21]),
     vec([82 22 10 70 49 40 34 48 80 71]),
     vec([91 17 62 75 47 11 7 72 35 55]),
     vec([71 90 75 64 94 15 12 67 20 50]),
     vec([70 93 77 29 58 93 68 57 7 52]),
     vec([87 63 26 6 82 27 56 48 36 95]),
     vec([36 15 41 78 76 84 30 76 36 8]),
     vec([88 81 13 82 54 13 29 40 78 75]),
     vec([88 54 64 32 52 6 54 82 6 26])]








     machines = [vec([3 4 6 5 1 8 9 10 2 7]),
     vec([5 8 2 9 1 4 3 6 10 7]),
     vec([10 7 5 4 2 1 9 3 8 6]),
     vec([2 3 8 6 9 5 4 7 10 1]),
     vec([7 2 4 1 3 9 5 8 10 6]),
     vec([8 6 9 3 5 7 4 2 10 1]),
     vec([7 2 5 6 3 4 8 9 10 1]),
     vec([1 6 9 10 4 7 5 8 3 2]),
     vec([6 3 4 7 5 8 9 10 2 1]),
     vec([10 5 7 8 1 3 9 6 4 2])]
    =#

    #=

     machines = [[2, 3, 1],
     [3, 1, 2],
     [3, 2, 1]]
    mach=machines
    machines = 0
     for i in 1:m
        a=hcat(vec(mach[i]),1:m)
        println(a)
        a=sortslices(a,dims=1)
        mach[i]=(a[:,2])
        end

    machines = mach

    println(machines)
    println(mach)
     =#

    model = Model(Gurobi.Optimizer)
    #model=Model(Cbc.Optimizer)      

    #set_optimizer_attribute(model, "logLevel", 1)
    #set_optimizer_attribute(model, "maxNodes", 100000)
    #set_optimizer_attribute(model, "seconds", 3600)
    @variable(model, c >= 0)
    @variable(model, y[1:n, 1:n, 1:m], Bin)
    @variable(model, x[1:n, 1:m] >= 0)






    M = sum(times[i, j] for i = 1:n for j = 1:m)
    #=
    for (j, i) in Iterators.product(1:n, 2:m)
     @constraint(model,x[j,machines[j][i]] - x[j,machines[j][i-1]] >= times[j][machines[j][i-1]])

     end
    =#

    for j = 1:n
        p = 1
        for i = 2:m
            @constraint(model, x[j, machines[j, i]] - x[j, machines[j, i-1]] >= times[j, p])
            p = p + 1
        end
    end

    times_old = zeros(n, m)
    for j = 1:n
        for i = 1:m
            val = times[j, i]
            times_old[j, machines[j, i]] = val
        end
    end


    for (j, k) in Iterators.product(1:n, 1:n)
        if k != j
            for i = 1:m
                @constraint(model, x[j, i] - x[k, i] + M * y[j, k, i] >= times_old[k, i])
                @constraint(
                    model,
                    -x[j, i] + x[k, i] - M * y[j, k, i] >= times_old[j, i] - M
                )
            end
        end
    end


    for j = 1:n
        @constraint(model, c - x[j, machines[j, m]] >= times[j, m])
    end

    @objective(model, Min, c)
    print(model)
    optimize!(model)


    #=                 
    for (j, i) in Iterators.product(1:n, 1:m)
       println(string("task ",j," starts on machine ",machines[j][i] ," at time", value(x[j,machines[j][i]]), " with time ", times[j][machines[j][i]]))
      end
    =#
    println("--------")

    for (j, i) in Iterators.product(1:n, 1:m)
        println(
            string(
                "task ",
                j,
                " starts on machine ",
                i,
                " at time",
                value(x[j, i]),
                " with time ",
                times_old[j, i],
            ),
        )
    end

    figure()

    for (j, i) in Iterators.product(1:n, 1:m)
        a = length(value(x[j, i]).+0.1:0.1:value(x[j, i])+times_old[j, i].-0.1)
        plot(value(x[j, i]).+0.1:0.1:value(x[j, i])+times_old[j, i].-0.1, i .* ones(a))
        text(value(x[j, i]) + times_old[j, i] / 2, i, string("job ", j))
    end




end

Run()
