function epiData = variedUptakeEpiSim(model, method, savename, rxnid, A, fredux, WTFlux, grData)

dlvl = 0.01
scaled = 0;

ngen = length(model.genes);
nrxn = length(model.rxns);
grRateKOTens = zeros(ngen, ngen, 11, 11);
grSingleDel = zeros(ngen,1);
grSingleMut = zeros(ngen,ngen,10);
grSingleMutrev = zeros(ngen,ngen,10);

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
epiData = zeros(lA,ngen,ngen);

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
    fluxtmp =  geometricFBA(mtmp);
    if numel(fluxtmp) < 1
      fluxtmp = optimizeCbModel(mtmp, 'max', 'one');
    end
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

if nargin < 8
  grData = zeros(lA,ngen,ngen);
  %put parfor back here later
  disp(lA)
  for i=1:lA
    mtmp = model;
    if A(i) < 0
      mtmp.lb(rxnid) = A(i);
    else
      mtmp.ub(rxnid) = A(i);
    end
    dtmp = doubleGeneMutationIsoAvg(mtmp, method, [fredux], [fredux], squeeze(WTFlux(i,:)), dlvl);
    grData(i,:,:) = dtmp;
    disp(i);
  end
  save(strcat(savename,'_gr.mat'),'grData');
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
  grSingleMut = singleGeneMutation(mtmp,method,squeeze(WTFlux(i,:)),dlvl);
  grSingleMutrev = zeros(ngen,ngen,10);
  for k = 1:ngen
    for j = 1:ngen
      grSingleMutrev(j,k,:) = grSingleMut(k,j,:);
    end
  end
  grRateKOTens(:,:,cellij,cellij) = grData(i,:,:);
  grRateKOTens(:,:,1:10,11)=grSingleMut;
  grRateKOTens(:,:,11,1:10)=grSingleMutrev;
  grRateKOTens = abs(grRateKOTens);
  grRateKOTens(:,:,11,11) = ones(ngen, ngen) * grRateWT(i);
  if scaled
    [intPos, intNeg ,eEffect,oInt] = findEpistaticDistributionsN11x11_scaled(mtmp,grRateKOTens,mtmp.genes,mtmp.genes,0.01);
  else
    [intPos, intNeg ,eEffect,oInt] = findEpistaticDistributionsN11x11(mtmp,grRateKOTens,mtmp.genes,mtmp.genes,0.01);
  end
  epiData(i,:,:) = squeeze(eEffect(:,:,cellij,cellij));
  dlmwrite(strcat(savename,num2str(A(i))), squeeze(grRateKOTens(:,:,cellij,cellij)));
  dlmwrite(strcat('epi',savename,num2str(A(i))), squeeze(epiData(i,:,:)));  
end

save(strcat(savename,'_epi.mat'),'epiData');
uptakeFlux=A;
save(strcat(savename,'_rxn',num2str(rxnid),'.mat'),'uptakeFlux');

