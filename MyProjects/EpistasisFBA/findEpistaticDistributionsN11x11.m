function [interactionsPos, interactionsNeg ,epistaticEffect,oppositeInt] = findEpistaticDistributionsN11x11(model,doubleDeletionFitness,geneList1,geneList2,minEffect)
%!!!!!!!!!!!! Read this !!!!!!!!!!!!!!!!!!!
  %Check for j =, depending on the type of gene matrix being used
  %Check the index being used as the wildtype growth value (usually, middle or last)
  %Set WT Fitness


%findEpistaticDistributions Finds synthetic lethal and/or synthetic sick
%interactions based on double deletion analysis data
%
  % [interactions,epistaticEffect,oppositeInt] = findEpistaticDistributions(model,doubleDeletionFitness,geneList1,geneList2,minEffect)
%
% model                 COBRA model
% doubleDeletionFitness A matrix of fitness (or growth rate) values for
%                       each of the double deletion strains. The diagonal
%                       of this matrix contains the single deletion fitness
%                       values.
% minEffect             Minimum fitness effect considered to be significant
%                       (opt, default 1e-2)
% 
% interactions          A sparse binary matrix indicating a SL or SS
%                       interaction between two genes in the model
% oppositeInt           A binary matrix indicating whether or not 
%                       there are epirstatic interactions of different sign
%                       occurring for some flux values.
% epistaticEffect       Magnitude of the epistatic interaction defined as
%                       f12-f1*f2
%                       %min(f1-f12,f2-f12) where f1 and f2 are the fitness
%                       values for the deletion strain of gene 1 and gene 2
%                       respectively and f12 is the fitness value for the
%                       double deletion strain of genes 1 and 2
% 
% The criteria for establishing a synthetic sick interaction are that the
% double deletion strain fitness must be at least minEffect lower than the
% fitness of either of the single deletion strains, i.e. 
%       f12 < f1-minEffect and f12 < f2-minEffect
%
% The additional criterion for establishing a synthetic lethal interaction
% is that the double deletion fitness value is smaller than minEffect (i.e.
% essentially zero)
%       f12 < minEffect
% 
% Note that the interactions matrix double counts all interactions
%
% Markus Herrgard 1/17/07

if (nargin < 3)
   geneList1 = model.genes;
end
if (nargin < 4)
  geneList2 = model.genes;
end 
if (nargin < 5)
    minEffect = 1e-2;
end

ZeroCut=4.5*10^-9;

%solMax = optimizeCbModel(constrainEnzFluxGlobal(model, 1),'max');
solMax = optimizeCbModel(model,'max');

grRateWT = solMax.f
disp(min(min(min(min(doubleDeletionFitness)))));
disp(max(max(max(max(doubleDeletionFitness)))));
doubleDeletionFitness = doubleDeletionFitness/grRateWT;
disp(min(min(min(min(doubleDeletionFitness)))));
disp(max(max(max(max(doubleDeletionFitness)))));

nGenes1 = length(geneList1);
nGenes2 = length(geneList2);
%singleDeletionFitness = diag(doubleDeletionFitness);

interactionsPos = zeros(nGenes1,nGenes2);
interactionsNeg = zeros(nGenes1,nGenes2);
oppositeInt = zeros(nGenes1,nGenes2);
epistaticEffect = zeros(nGenes1,nGenes2,11,11);

for i = 1:nGenes1
    %[isInModel,geneID1] = ismember(geneList1{i},model.genes);
    for j = i+1:nGenes2
       %[isInModel,geneID2] = ismember(geneList2{j},model.genes);
       posC = 0;
       negC = 0;
	   for mval1 = 1:11
          for mval2 = 1:11
			  fitness1 = doubleDeletionFitness(i,j,mval1,11);
              fitness2 = doubleDeletionFitness(i,j,11,mval2);
              fitness12 = doubleDeletionFitness(i,j,mval1,mval2);
              if fitness1 < ZeroCut
			    fitness1 = 0;
			  end
              if fitness2 < ZeroCut
                fitness2 = 0;
			  end
              if fitness12 < ZeroCut
                fitness12=0;
			  end			  
			  %fitness1 = doubleDeletionFitness(i,j,mval1,11);
              %fitness2 = doubleDeletionFitness(i,j,11,mval2);
              %fitness12 = doubleDeletionFitness(i,j,mval1,mval2);
              epistaticEffect(i,j,mval1,mval2) = fitness12-fitness1*fitness2;
              epistaticEffect(j,i,mval2,mval1) = epistaticEffect(i,j,mval1,mval2);
              if epistaticEffect(i,j,mval1,mval2) >= minEffect
                 interactionsPos(i,j) = 1;
              elseif epistaticEffect(i,j,mval1,mval2) <= -1*minEffect
	             interactionsNeg(i,j) = -1;
              end
              if epistaticEffect(i,j,mval1,mval2) > 0 & abs(epistaticEffect(i,j,mval1,mval2)) >= minEffect
                 posC = posC + 1;
              end
              if epistaticEffect(i,j,mval1,mval2) < 0 & abs(epistaticEffect(i,j,mval1,mval2)) >= minEffect
                 negC = negC + 1;
              end
           end
        end
        if posC > 0 & negC > 0
           oppositeInt(i,j) = 1;
        end
        %make histograms or heatmaps here?
    end
end
