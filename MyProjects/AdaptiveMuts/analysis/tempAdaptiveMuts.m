function [genesTadapted, mutLoss] = tempAdaptiveMuts(AdaptedFitness, RelativeFitness) 
%Requires output from makeSDTrajectory.m
%Look at WT growth rates that are at least 1% of maximum:
[tpoints, ngen] = size(RelativeFitness);

adapt_start = -1;
genesTadapted = zeros(1,ngen);
mutLoss = zeros(1,ngen);

for i = 1:tpoints
  if AdaptedFitness(i,1)/AdaptedFitness(tpoints,1) >= 0.01
    adapt_start = i;
    break;
  end
end

for i = 1:ngen
  if max(RelativeFitness(adapt_start:tpoints,i)) >= 1.001
    genesTadapted(i) = 1;  
  end
end

for i = 1:ngen
  if AdaptedFitness(tpoints,i+1) < AdaptedFitness(tpoints,1) * 0.9999
    if genesTadapted(i) > 0
      mutLoss(i) = 1;
    end
  end
end