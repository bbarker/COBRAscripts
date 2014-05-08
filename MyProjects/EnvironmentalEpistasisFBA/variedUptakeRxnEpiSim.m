function epiData = variedUptakeRxnEpiSim(model, method, savename, rxnid, A, fredux, WTFlux, grData)
%!!!Consider saving other phenotypes after this works
%!!!This may take some constructive saving to avoid memory issues:
%perhaps, create a directory or file for each cell and save the simulation result there
%dlvl = 0.01

geneRxns = any(model.rxnGeneMat,2);
ngrxn = sum(geneRxns);
ngen = length(model.genes);
nrxn = length(model.rxns);
grRateKOTens = zeros(ngrxn, ngrxn, 11, 11);
grSingleDel = zeros(ngrxn,1);
grSingleMut = zeros(ngrxn,ngrxn,10);
grSingleMutrev = zeros(ngrxn,ngrxn,10);

%make params:
%A = -1*[0.1 0.5 [1:20]]; %glucose uptake lvl, check lowest is "feasible"
if nargin < 5
  A = -9 - [2:9]/10;
end

if nargin < 6
  fredux = 0;
end

cellij = floor(10*fredux+1.0001);
lA = length(A);
grRateWT = zeros(lA,1);
epiData = zeros(lA,ngrxn,ngrxn);

if nargin < 7
  WTFlux = zeros(lA,nrxn);
  parfor i=1:lA
    mtmp = model;
    if A(i) < 0
      mtmp.lb(rxnid) = A(i);
    else
      mtmp.ub(rxnid) = A(i);
    end
    solt = optimizeCbModel(mtmp,'max');
    if abs(solt.f) < 1e-7
      error(strcat('No solution found for problem ',num2str(i)))
    end 
    grRateWT(i) = solt.f;
    mtmp.lb(find(model.c)) = solt.f; %Only works for single-valued objective functions
    disp(strcat('Starting geometric WT sim ',num2str(i)));
    fluxtmp =  geometricFBA2(mtmp);
    if norm(fluxtmp,1) < 1e-7
      sol1n = optimizeCbModel(mtmp,'max','one'); 
      WTFlux(i,:) = (sol1n.x)';
    else
      WTFlux(i,:) = fluxtmp;
    end
    disp(strcat('Finished geometric WT sim ',num2str(i)));
  end
  potentialFail = abs(WTFlux(:,[536,1577])) < 1e-7;
  pFsum = sum(potentialFail(:,1) & potentialFail(:,2));
  save(strcat(savename,'_Flux.mat'),'WTFlux');
  if pFsum > 0
    error(strcat(num2str(pFsum)),' unlikely solutions encountered from geoFBA')
  end
end


%just in case:
changeCobraSolver('gurobi','LP');

if nargin < 8
  grData = zeros(lA,ngrxn,ngrxn);
  parfor i=1:lA
    mtmp = model;
    if A(i) < 0
      mtmp.lb(rxnid) = A(i);
    else
      mtmp.ub(rxnid) = A(i);
    end
    dtmp = doubleRxnMutation(mtmp, method, [fredux], [fredux], squeeze(WTFlux(i,:)));
    grData(i,:,:) = dtmp;
 disp(i);
  end
  save(strcat(savename,'_Rxngr.mat'),'grData');
end

%We'll want to save the grData mat, WT fluxes, and single Mutant fit vectors, 
%but also compute epistasis and save the array of single cells by substituting in to a "sparse" grTens mat 
%so we don't need to save such large matrices for no reason

disp('Starting Epistasis  Loop');
for i=1:lA
  mtmp = model;
    if A(i) < 0
      mtmp.lb(rxnid) = A(i);
    else
      mtmp.ub(rxnid) = A(i);
    end
  %Could optimize the single-mutation code below ... but premature optimization and all that
  disp(strcat('Entering single mutant calc ', num2str(i)));
  grSingleMut = singleRxnMutation(mtmp,method,squeeze(WTFlux(i,:)));
  grSingleMutrev = zeros(ngrxn,ngrxn,10);
  for k = 1:ngrxn
    for j = 1:ngrxn
      grSingleMutrev(j,k,:) = grSingleMut(k,j,:);
    end
  end
  grRateKOTens(:,:,cellij,cellij) = grData(i,:,:);
  grRateKOTens(:,:,1:10,11)=grSingleMut;
  grRateKOTens(:,:,11,1:10)=grSingleMutrev;
  grRateKOTens = abs(grRateKOTens);
  grRateKOTens(:,:,11,11) = ones(ngrxn, ngrxn) * grRateWT(i);
  [intPos, intNeg ,eEffect,oInt] = findEpistaticDistributionsN11x11(mtmp,grRateKOTens,mtmp.rxns(geneRxns),mtmp.rxns(geneRxns),0.01);
  epiData(i,:,:) = squeeze(eEffect(:,:,cellij,cellij));
  dlmwrite(strcat(savename,num2str(A(i))), squeeze(grRateKOTens(:,:,cellij,cellij)));
  dlmwrite(strcat('epiRxn',savename,num2str(A(i))), squeeze(epiData(i,:,:)));  
end

save(strcat(savename,'_epiRxn.mat'),'epiData');
uptakeFlux=A;
save(strcat(savename,'_rxn',num2str(rxnid),'.mat'),'uptakeFlux');

