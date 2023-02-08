using DelimitedFiles
using Random


function run()




PT = [[1, 2, 2],
         [2, 1, 2],
         [1, 2, 1]]

Ma = [[3, 1, 2],
            [2, 3, 1],
            [3, 2, 1]]


            j_num = 3

            ma_num = 3

population_size = 30
crossover_rate = 0.9
mutation_rate = 0.1
Num_Iteration = 1000


population_list = zeros(population_size, j_num * ma_num)

for m in 1:population_size
    population_list[m,:] = shuffle(collect(1:j_num * ma_num))
    for j in 1:j_num * ma_num
        for k in 0:(j_num - 1)
            if population_list[m, j] > k * ma_num && population_list[m, j] <= (k + 1) * ma_num
                population_list[m, j] = k + 1
            end
        end
    end
end


Makespan_best = 9999999999; 
Scheduling_best = zeros(1, j_num*ma_num);  


for i in 1:Num_Iteration


    # Crossover
    population_list_tmp = copy(population_list)
    S = shuffle(collect(1:population_size))
    for m in 1:(population_size÷2)
        crossover_prob = rand()
        if (crossover_rate >= crossover_prob) 
            parent_1 = population_list[S[2 * m - 1], :]
            parent_2 = population_list[S[2 * m], :]
            child_1 = copy(parent_1)
            child_2 = copy(parent_2)
            cutpoint = sort(randperm(j_num * ma_num)[1:2])
            for k in cutpoint[1]:cutpoint[2]
                child_1[k] = parent_2[k]
                child_2[k] = parent_1[k]
            end
            population_list[S[2 * m - 1], :] = child_1
            population_list[S[2 * m], :] = child_2
        end
    end


    # Mutation
for m = 1:population_size
    for j = 1:j_num * ma_num
        mutation_prob = rand()
        if mutation_rate >= mutation_prob
            ran_num = rand()
            for k = 0:(j_num - 1)
                if ran_num > k * (1 / j_num) && ran_num <= (k + 1) * (1 / j_num)
                    population_list[m, j] = k + 1
                end
            end
        end
    end
end

# Repairment
for m = 1:population_size
    repair_or_not = zeros(1, j_num)
    for j = 1:j_num * ma_num
        for k = 1:j_num
            if population_list[m, j] == k
                repair_or_not[k] += 1
            end
        end
    end

    for k = 1:j_num
        if repair_or_not[k] > ma_num


            r_ran_num = shuffle(collect(1:repair_or_not[k])) 
            println(r_ran_num[1:Int64(repair_or_not[Int64(k)]-ma_num)])
            r_ran_num = sort(r_ran_num[1:Int64(repair_or_not[Int64(k)]-ma_num)])     




            #r_ran_num = sort(randperm(repair_or_not[k])[1:(repair_or_not[k] - ma_num)])
            appeartime = 0
            for j = 1:j_num * ma_num
                if population_list[m, j] == k
                    appeartime += 1
                    for n = 1:(repair_or_not[k] - ma_num)
                        if appeartime == r_ran_num[Int64(n)]
                            population_list[m, j] = 0
                        end
                    end
                end
            end
        end
    end

    for k = 1:j_num
        if repair_or_not[k] < ma_num
            zeroappeartime = 0
            for j = 1:j_num * ma_num
                if population_list[m, j] == 0
                    zeroappeartime += 1
                end
            end

            r_ran_num = shuffle(collect(1:zeroappeartime)) 
            r_ran_num = sort(r_ran_num[1:Int64(ma_num-repair_or_not[Int64(k)])])


            #r_ran_num = sort(randperm(zeroappeartime)[1:ma_num - (repair_or_not[k])])
            appeartime = 0
            for j = 1:j_num * ma_num
                if population_list[m, j] == 0
                    appeartime += 1
                    for n = 1:(ma_num - repair_or_not[k])
                        if appeartime == r_ran_num[Int64(n)]
                            population_list[m, j] = k
                        end
                    end
                end
            end
        end
    end
end



population_list = [population_list_tmp; population_list]

fitness = zeros(population_size*2, 1+j_num*ma_num)
Gen = zeros(population_size*2,j_num*ma_num)
Gen_m = zeros(population_size*2,j_num*ma_num)
Gen_t = zeros(population_size*2,j_num*ma_num)  
MachineJob = zeros(population_size*2,j_num*ma_num)
MachineTimeBegin = zeros(population_size*2,j_num*ma_num) 
MachineTimeEnd = zeros(population_size*2,j_num*ma_num) 

for m in 1:population_size*2 
    Gen[m,1:j_num*ma_num] = population_list[m,1:j_num*ma_num]
    Ma2 = Ma
    PT2 = PT
    for j in 1:j_num*ma_num
        for k in 1:ma_num
            if Ma2[Int64(Gen[m,j]),k] != 0
                Gen_m[m,j] = Ma2[Int64(Gen[m,j]),k]
                Gen_t[m,j] = PT2[Int64(Gen[m,j]),k]
                Ma2[Int64(Gen[m,j]),k] = 0
                break
            end
        end
    end

    t = 1 
    for k in 1:ma_num
        for j in 1:j_num*ma_num
            if Gen_m[m,j] == k
                MachineJob[m,t] = Gen[m,j] 
                t = t + 1
            end
        end
    end
end




for m in 1:population_size*2
    j_count = zeros(1,j_num)
    m_count = zeros(1,ma_num)
    t_count = zeros(1,ma_num)

    for j in 1:j_num*ma_num
        revise = 0
        t_count[Gen_m[m,j]] += 1
        if t_count[Gen_m[m,j]] <= 2
            j_count[Gen(m, j)] += Gen_t(m, j)
            m_count[Gen_m(m, j)] += Gen_t(m, j)
                if m_count[Gen_m(m, j)] < j_count[Gen(m, j)]
                    m_count[Gen_m(m, j)] = j_count[Gen(m, j)]
                elseif m_count[Gen_m(m, j)] > j_count[Gen(m, j)]
                    j_count[Gen(m, j)] = m_count[Gen_m(m, j)]
                end
            MachineTimeBegin[m, Gen_m(m, j) * j_num - j_num + t_count(Gen_m(m, j))] = m_count[Gen_m(m, j)] - Gen_t(m, j)
            MachineTimeEnd[m, Gen_m(m, j) * j_num - j_num + t_count(Gen_m(m, j))] = m_count[Gen_m(m, j)]





        elseif t_count[Gen_m[m,j]] >= 3
            temBeginEnd = zeros(2,t_count[Gen_m[m,j]]);

            for k = 1:(t_count[Gen_m[m,j]]-1)
                temBeginEnd[1, :] = MachineTimeBegin[m, Gen_m(m, j) * j_num - j_num + 1 : Gen_m(m, j) * j_num - j_num + t_count(Gen_m(m, j)) - 1]
                temBeginEnd[2, :] = MachineTimeEnd[m, Gen_m(m, j) * j_num - j_num + 1 : Gen_m(m, j) * j_num - j_num + t_count(Gen_m(m, j)) - 1]
                temBeginEnd = sortrows(temBeginEnd', dims = 1)
                temBeginEnd = transpose(temBeginEnd)

                if j_count[Gen[m,j]] <= temBeginEnd[1,k+1]
                    if j_count[Gen[m,j]] >= temBeginEnd[2,k]
                        if temBeginEnd[1,k+1] - j_count[Gen[m,j]] >= Gen_t[m,j]
                            MachineTimeBegin[m,Gen_m[m,j]*j_num-j_num+t_count[Gen_m[m,j]]] = j_count[Gen[m,j]]
                            MachineTimeEnd[m,Gen_m[m,j]*j_num-j_num+t_count[Gen_m[m,j]]] = j_count[Gen[m,j]] + Gen_t[m,j];
                            break
                        end
                    elseif j_count[Gen[m,j]] < temBeginEnd[2,k]
                        if temBeginEnd[1,k+1] - temBeginEnd[2,k] >= Gen_t[m,j]
                            revise = 1;
                            MachineTimeBegin[m,Gen_m[m,j]*j_num-j_num+t_count[Gen_m[m,j]]] = MachineTimeEnd[m,Gen_m[m,j]*j_num-j_num+k]
                            MachineTimeEnd[m,Gen_m[m,j]*j_num-j_num+t_count[Gen_m[m,j]]] = MachineTimeEnd[m,Gen_m[m,j]*j_num-j_num+k] + Gen_t[m,j]
                            break
                        end
                    end
                end                
            end


            if revise == 0
                j_count[Gen[m,j]] = j_count[Gen[m,j]] + Gen_t[m,j]
                m_count[Gen_m[m,j]] = m_count[Gen_m[m,j]] + Gen_t[m,j]
                if m_count[Gen_m[m,j]] < j_count[Gen[m,j]]
                    m_count[Gen_m[m,j]] = j_count[Gen[m,j]]
                elseif m_count[Gen_m[m,j]] > j_count[Gen[m,j]]
                    j_count[Gen[m,j]] = m_count[Gen_m[m,j]]
                end
                MachineTimeBegin[m,Gen_m[m,j]*j_num-j_num+t_count[Gen_m[m,j]]] = m_count[Gen_m[m,j]] - Gen_t[m,j]
                MachineTimeEnd[m,Gen_m[m,j]*j_num-j_num+t_count[Gen_m[m,j]]] = m_count[Gen_m[m,j]]               
            elseif revise == 1
                j_count[Gen[m,j]] = MachineTimeEnd[m,Gen_m[m,j]*j_num-j_num+t_count[Gen_m[m,j]]]
            end




        end
        fitness[m,1] = max[j_count]
        fitness[m,2:1+j_num*ma_num] = Gen[m,1:j_num*ma_num]

    end
                        
    Totalfitness = 0
    for m = 1:population_size
        Totalfitness = Totalfitness + 1/fitness[m,1]
    end
    
    pk = zeros(population_size,1)
    for m = 1:population_size
        pk[m] = [1/fitness[m,1]] / Totalfitness
    end
    
    qk = zeros(population_size,1)
    for m = 1:population_size 
        cumulative = 0
        for j = 1:m
            cumulative = cumulative + pk[j]
        end
        qk[m] = cumulative
    end
    
    last_population_list = population_list
    population_list = zeros(population_size, j_num*ma_num)
    
    selection_rand = zeros(1,population_size)
    for m = 1:population_size
        selection_rand[m] = rand()  ############################################ fix this
    end


    for m = 1:population_size
        if selection_rand[m] <= qk[1]
            population_list[m, 1:j_num*ma_num] = last_population_list[m, 1:j_num*ma_num]
        else
            for j = 1:(population_size-1)
                if selection_rand[m] > qk[j] && selection_rand[m] <= qk[j+1]
                    population_list[m, 1:j_num*ma_num] = last_population_list[(j+1), 1:j_num*ma_num]
                    break
                end
            end
        end
    end


end

end

end


run()
