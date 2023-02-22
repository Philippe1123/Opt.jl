
using XLSX





function run()


    nameOfFile = "/home/philippe/.julia/dev/Opt/src/JSP_dataset_small.xlsx"
    xf = XLSX.readxlsx(nameOfFile)
    n = xf["Parameters!B2"]
    m = xf["Parameters!B3"]



    times = xf[string("ProcessingTime!B2:", xf[string("Sheet4!A", m + 1)], n + 1)]

    machines = xf[string("MachineSequence!B2:", xf[string("Sheet4!A", m + 1)], n + 1)]



    mat = Matrix(transpose(reshape(2:n*m+1, n, m)))
    GraphMat = zeros(n * m + 2, n * m + 2)
    TimesMat = zeros(n * m + 2, n * m + 2)
    TimesDict = Dict()

    for j = 1:n
        for i = 1:m
            TimesDict[mat[i,j]] =  times[i,j]
        end
    end

    for j = 1:n
        for i = 1:m-1

            GraphMat[mat[j, i], mat[j, i+1]] = 1
            TimesMat[mat[j, i], mat[j, i+1]] = times[j,i]
            println(GraphMat[mat[j, i], mat[j, i+1]])

        end
    end

    for j = 1:n

        GraphMat[mat[j, m], n*m+2] = 1
        TimesMat[mat[j, m], n*m+2] = times[j,end]
        #println(GraphMat[mat[j,m],n*m+2])


    end


    for j = 1:n

        GraphMat[1, mat[j, 1]] = 1
        #println(1,GraphMat[mat[j,m]])


    end

    for j = 1:n-1
        for i = 1:m
            out = machines[j, i]

            for p = j+1:n
                for k = 1:m
                    println([j, i, p, k])
                    if (out == machines[p, k])
                        GraphMat[mat[j, i], mat[p, k]] = 2
                        GraphMat[mat[p, k], mat[j, i]] = 2
                    end
                end
            end
            println("---")


        end
    end



    @show GraphMat

    Bidir(GraphMat, TimesMat,TimesDict, mat, 3)




end




function Bidir(Graph::Matrix, TimesMat::Matrix,TimesDict::Dict ,Jobs::Matrix, njobs::Int64)
    TimesMatInv = deepcopy(TimesMat)
    #initialisation
    Sol = deepcopy(Graph)
    L = Int64[]
    append!(L, 1)
    R = Int64[]
    append!(R, size(Sol, 2))
    r = Dict()
    t = Dict()
    S = Int64[]
    T = Int64[]
    pos = findall(x -> x == 1, Graph[1, :])
    ln = length(pos)
    for op = 1:ln
        append!(S, pos[op])
        r[pos[op]] = 0 
    end

    pos = findall(x -> x == 1, Graph[:, end])
    ln = length(pos)
    for op = 1:ln
        append!(T, pos[op])
        t[pos[op]] = 0
    end
    SJ = Dict()
    PJ = Dict()
    for j = 1:size(Jobs, 1)
        for i = 1:size(Jobs, 2)-1
            SJ[Jobs[j, i]] = Jobs[j, i+1]
        end

        for i = size(Jobs, 2):-1:2
            PJ[Jobs[j, i]] = Jobs[j, i-1]
        end
    end

    #njobs = 11# change

    while (length(L) + length(R)) != size(Graph, 2) 
        n_in_S = length(S)
        if n_in_S != 0
            elem = 1
        end
        operation_under_consideration = S[elem]

        A = Sol[operation_under_consideration, :]
        pos = findall(x -> x == 2, A)

        ln = length(pos)

        for op = 1:ln
            Sol[operation_under_consideration, pos[op]] = 1
            TimesMat[pos[op], operation_under_consideration]=TimesMat[pos[op], operation_under_consideration]+TimesDict[operation_under_consideration]
            Sol[pos[op], operation_under_consideration] = 0
        end
        append!(L, operation_under_consideration)
        remove!(T, operation_under_consideration)
        remove!(S, operation_under_consideration)


        pos = findall(x -> x == SJ[operation_under_consideration], R)
        if length(pos) == 0
            append!(S, SJ[operation_under_consideration])
            r[SJ[operation_under_consideration]] = 0 # add new operation to r
            
        end


        # Update r
        for el = 1 :length(S)
            int = Sol[:,S[el]]
            pos = findall(x -> x == 1, int)
            for p = 1:length(pos)
            r[S[el]] = r[S[el]] + TimesMat[S[el],pos[p]] # machine precedence constraint
            TimesMat[S[el],pos[p]] = 0
            r[S[el]] = r[S[el]] + TimesMat[pos[p],S[el]] # job precedence constraint
            TimesMat[pos[p],S[el]] = 0
            end
        end

        


        if (length(L) + length(R)) != size(Graph, 2)
            n_in_S = length(S)
            if n_in_S != 0
                elem = 1
            end
            operation_under_consideration = T[elem]
            A = Sol[operation_under_consideration, :]
            pos = findall(x -> x == 2, A)

            ln = length(pos)

            for op = 1:ln
                Sol[operation_under_consideration, pos[op]] = 0
                #TimesMatInv[pos[op], operation_under_consideration]=TimesMatInv[pos[op], operation_under_consideration]+TimesDict[operation_under_consideration]#check
                TimesMatInv[operation_under_consideration, pos[op]]=TimesMatInv[operation_under_consideration, pos[op]]+TimesDict[operation_under_consideration]#check

                Sol[pos[op], operation_under_consideration] = 1

                
            end

            append!(R, operation_under_consideration)
            remove!(T, operation_under_consideration)
            remove!(S, operation_under_consideration)

            pos = findall(x -> x == PJ[operation_under_consideration], R)
            if length(pos) == 0
                append!(T, PJ[operation_under_consideration])
                t[PJ[operation_under_consideration]] = 0 # add new operation to t #check
            end
            # should be an update of t now
            #check whole block
            for el = 1 :length(T)
                int = Sol[T[el],:]
                pos = findall(x -> x == 1, int)
                for p = 1:length(pos)
                t[T[el]] = t[T[el]] + TimesMatInv[T[el],pos[p]] # machine precedence constraint
                TimesMatInv[T[el],pos[p]] = 0
                t[T[el]] = t[T[el]] + TimesMatInv[pos[p],T[el]] # job precedence constraint
                TimesMatInv[pos[p],T[el]] = 0
                end
            end


        end


    end
    println(L)
    println(R)
    println(r)
    println(t)



    #=
            println(a)
            a = a + 1

        end

        for j = 1:njobs
            for i = 1:noperations


                A = Sol[Jobs[j, i], :]

                pos = findall(x -> x == 2, A)

                ln = length(pos)
                println(Jobs[j, i])


                for op = 1:ln

                    Sol[S[j][i], pos[op]] = 1
                    Sol[pos[op], S[j][i]] = 0
                    #println(pos[op])
                    #Sol[j, pos[i_inv]] = 1
                    #Sol[pos[i_inv], j] = 0
                    #i_inv = i_inv - 1
                end
                append!(L, S[j][i])
                remove!(T[j], S[j][i])
                remove!(S[j], S[j][i])

                if (length(L) + length(R) != size(Sol, 2))

                    println("Cont")
                end
            end

        end
        println(L)
        println(T)
        println(S)
        println("ok")
        @show Sol

    =#

end


function remove!(a::Array{}, item::Int64)
    deleteat!(a, findall(x -> x == item, a))
end




run()
