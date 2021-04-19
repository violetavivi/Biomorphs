% clear all;

populationSize = 9;
numberOfGenes = 12;
% n,rR, RL,rSigma,phiR, phiL,phiSigma,x0,y0,colorTreeR,colorTreeG,colorTreeB
lengthGenes = 20;
crossoverProbability = 0.8;
mutationProbability = 0.025;
tournamentSelectionParameter = 0.8;
variableRangeMatrix = [4 7; %n
    0.2 0.9; % rR
    0.2 0.9; % rL
    0   0.1; % rSigma
    -pi/2 pi/2;  % phiR
    -pi/2 pi/2;  % phiL
    0   0.1; % phiSig
    -0.1 0.1;   % x1
    0.2 1;   % y1
    0 1;   %colorTree
    0 1;
    0 1];

numberOfGenerations = 100;
x0 = 0;
y0 =  0;
population = InitializePopulation(populationSize,lengthGenes, numberOfGenes);
%%
for iGeneration = 1:numberOfGenerations %generation loop
    
    maximunFitness = 0.0;
    xBest = zeros(1,2);
    
    %--------------------- EVALUATION (aka compute fitness)
    figure
    
    for ii = 1:populationSize
        ss = subplot(3,3,ii);
        values = DrawTree(reshape(population(ii,:,:),[lengthGenes, numberOfGenes])',x0,y0,ss,variableRangeMatrix);
        set(gca,'tag',num2str(ii))
    end
    sgtitle(sprintf('generation %d',iGeneration))
    bestIndividualIndex  = str2double(clicksubplot_withOutput);
    sgtitle(sprintf('generation %d, evoling...',iGeneration))
    distanceFromChosen = zeros(populationSize,1);
    for ii = 1:populationSize
        distanceFromChosen(ii) = norm(values(ii)-values(bestIndividualIndex));
    end
    distanceFromChosenNorm = distanceFromChosen/sum(distanceFromChosen);
    fitness = 1 - distanceFromChosenNorm;
    %%
    %----------- REPRODUCTION
    tempPopulation = population;
    
    for i = 1:2:populationSize
        i1 = TournamentSelect(fitness,tournamentSelectionParameter);
        i2 = TournamentSelect(fitness,tournamentSelectionParameter);
        
        for gg = 1:numberOfGenes
            chromosome1 = population(i1,:,gg);
            chromosome2 = population(i2,:,gg);
            
            r = rand;
            if (r < crossoverProbability)
                newChromosomePair = Cross(chromosome1,chromosome2);
                tempPopulation(i,:,gg) = newChromosomePair(1,:);
                tempPopulation(i+1,:,gg) = newChromosomePair(2,:);
            else
                tempPopulation(i,:,gg) = chromosome1;
                tempPopulation(i+1,:,gg) = chromosome2;
            end
        end
    end
    %%
    %----------- MUTATION
    for i = 1:populationSize
        for gg = 1:numberOfGenes
            originalChromosome = tempPopulation(i,:,gg);
            mutatedChromosome = Mutate(originalChromosome,mutationProbability);
            tempPopulation(i,:,gg) = mutatedChromosome;
        end
    end
    
    tempPopulation(1,:,:) = population(bestIndividualIndex,:,:);
    population = tempPopulation;


end % loop over generations

