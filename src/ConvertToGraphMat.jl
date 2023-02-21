
using XLSX





function run()


    nameOfFile = "JSP_dataset_small.xlsx"
    xf = XLSX.readxlsx(nameOfFile)
    n = xf["Parameters!B2"]
    m = xf["Parameters!B3"]



    times = xf[string("ProcessingTime!B2:", xf[string("Sheet4!A", m + 1)], n + 1)]

    machines = xf[string("MachineSequence!B2:", xf[string("Sheet4!A", m + 1)], n + 1)]



    mat = Matrix(transpose(reshape(2:n*m+1, n, m)))
    GraphMat = zeros(n * m + 2, n * m + 2)

    for j = 1:n
        for i = 1:m-1

            GraphMat[mat[j, i], mat[j, i+1]] = 1
            println(GraphMat[mat[j, i], mat[j, i+1]])

        end
    end

    for j = 1:n

        GraphMat[mat[j, m], n*m+2] = 1
        #println(GraphMat[mat[j,m],n*m+2])


    end


    for j = 1:n

        GraphMat[1, mat[j, m]] = 1
        #println(1,GraphMat[mat[j,m]])


    end



    #=
        for j = 1 : n - 1
            for i = 1 :m
                out = machines[j,i]
                for k = 1 :m
                if (out == machines[j+1,k])
                    println([mat[j,i],mat[j+1,k]])
                    println([mat[j+1,k],mat[j,i]])
                    GraphMat[mat[j,i],mat[j+1,k]] = 1  
                    GraphMat[mat[j+1,k],mat[j,i]] = 1  
                    #println("hit")

                end
            end
            println("---")


            end
        end 
    =#

    for j = 1:n-1
        for i = 1:m
            out = machines[j, i]

            for p = j+1:n
                for k = 1:m
                    println([j, i, p, k])
                    if (out == machines[p, k])
                        #println([mat[j,i],mat[p,k]])
                        #println([mat[p,k],mat[j,i]])
                        GraphMat[mat[j, i], mat[p, k]] = 2
                        GraphMat[mat[p, k], mat[j, i]] = 2
                        #println("hit")

                    end
                end
            end
            println("---")


        end
    end



    @show GraphMat

    Bidir(GraphMat, mat, 3)




end




function Bidir(Graph::Matrix, Jobs::Matrix, njobs::Int64)

    Sol = deepcopy(Graph)

    L = Int64[]
    append!(L, 1)

    R = Int64[]

    append!(R, size(Sol, 2))
    i_inv = njobs
    noperations = 3

    #S = Jobs[j, i]
    #T = Jobs[j, i_inv]

    S = Array{Int64}[]
    T = Array{Int64}[]

    PJ = Dict()
    for i = 1:njobs
        intermDict = Dict()
        for j = 1 : 3
            intermDict[]


        end



    end



    for j = 1 :njobs
        push!(S,Jobs[j,:])
        push!(T,reverse(Jobs[j,:]))
    end


    println(S)
    println(T)



    #njobs = 11# change

    for j = 1:njobs
        for i = 1: noperations


            A = Sol[Jobs[j,i], :]

            pos = findall(x -> x == 2, A)

            ln = length(pos)
            println(Jobs[j,i])


            for op = 1:ln

                Sol[S[j][i], pos[op]] = 1
                Sol[pos[op], S[j][i]] = 0
                #println(pos[op])
                #Sol[j, pos[i_inv]] = 1
                #Sol[pos[i_inv], j] = 0
                #i_inv = i_inv - 1
            end
            append!(L, S[j][i])
            remove!(T[j],S[j][i])
            remove!(S[j],S[j][i])

            if (length(L) + length(R) != size(Sol,2))

                println("Cont")
            end
        end

    end
    println(L)
    println(T)
    println(S)
    println("ok")
    @show Sol



end


function remove!(a::Array{}, item::Int64)
    deleteat!(a, findall(x -> x == item, a))
end




run()
