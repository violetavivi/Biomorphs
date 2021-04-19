function [x, xNorm] = DecodeChromosome(singleChromosome,singleVariableRange)
%%
nGenes = size(singleChromosome,2);

xNorm = 0.0;
for j = 1:nGenes
    xNorm = xNorm + singleChromosome(j)*2^(-j);
end
x = singleVariableRange(1) + (singleVariableRange(2)- singleVariableRange(1))*xNorm(1);


end
