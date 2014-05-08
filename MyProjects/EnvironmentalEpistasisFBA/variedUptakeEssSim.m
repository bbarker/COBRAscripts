function [grRateKOTens, essSums] = variedUptakeEssSim(model, method, savename, rxnid, A)
%!!!Consider saving other phenotypes after this works
%!!!This may take some constructive saving to avoid memory issues:
%perhaps, create a directory or file for each cell and save the simulation result there
dlvl = 0.01


ngen = length(model.genes);
nrxn = length(model.rxns);

%make params:
%A = -1*[0.1 0.5 [1:20]]; %glucose uptake lvl, check lowest is "feasible"
if nargin < 5
  A = -9 - [2:9]/10;
end

lA = length(A);
grRateKOTens = zeros(lA, ngen);
essSums = zeros(lA,1);

parfor i=1:lA
  mtmp = model;
  mtmp.lb(rxnid) = A(i);
  [grRatio,grSingleDel,grRateWT,hasEffect,delRxns,fluxSolution] = singleGeneDeletion(mtmp,method);
  %Could optimize the single-mutation code below ... but premature optimization and all that
  grRateKOTens(i,:) = grRatio;
  essSums(i) = sum(grRatio < (1e-5));
end

save(strcat(savename,'_SingleKOGradient.mat'),'grRateKOTens');
dlmwrite(strcat(savename,'_SingleKOGradient.csv'),grRateKOTens);
dlmwrite(strcat(savename,'_EssentialSum.csv'),essSums);
uptakeFlux=A;
save(strcat(savename,'_rxn',num2str(rxnid),'.mat'),'uptakeFlux');

