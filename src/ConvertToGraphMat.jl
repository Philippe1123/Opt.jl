
using XLSX





function run()


    nameOfFile = "JSP_dataset_small.xlsx"
    xf = XLSX.readxlsx(nameOfFile)
    n = xf["Parameters!B2"]
    m = xf["Parameters!B3"]
    
    
    
    times = xf[string("ProcessingTime!B2:", xf[string("Sheet4!A", m + 1)], n + 1)]
    
    machines = xf[string("MachineSequence!B2:", xf[string("Sheet4!A", m + 1)], n + 1)]



    mat = transpose(reshape(2:n*m+1,n,m))
    GraphMat = zeros(n*m+2,n*m+2)

    for j = 1 : n
        for i = 1 :m-1

            GraphMat[mat[j,i],mat[j,i+1]]=1
            println(GraphMat[mat[j,i],mat[j,i+1]])

        end
    end

    for j = 1 : n

            GraphMat[mat[j,m],n*m+2]=1
            println(GraphMat[mat[j,m],n*m+2])

        
    end


    for j = 1 : n

        GraphMat[mat[j,m],1]=1
        println(GraphMat[mat[j,m],1])

    
    end 




    for j = 1 : n - 1
        for i = 1 :m
            out = machines[j,i]
            for k = 1 :m
            if (out == machines[j+1,k])
                #println(mat[j+1,i])
                #println(mat[j+1,k])
                GraphMat[mat[j,i],mat[j+1,k]] = 1  
                GraphMat[mat[j+1,k],mat[j,i]] = 1  
                #println("hit")

            end
        end



        end
    end 



    @show GraphMat;




end




run()
