function population = InitializePopulation(populationSize,lenGenes, nGenes)

population = zeros(populationSize, lenGenes,nGenes);
for i = 1:populationSize
    for j = 1:lenGenes
        for z= 1:nGenes
            s = rand;
            if(s < 0.5)
                population(i,j,z) = 0;
            else
                population(i,j,z) = 1;
            end
        end
        
    end
end