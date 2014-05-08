function WTFlux = variedUptakeQuickWT(model, method, savename, rxnid, A, fredux)
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

if nargin < 6
  fredux = 0;
end

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
    mtmp.lb(nrxn) = solt.f; %May be different in some models or objectives
    sol1n = optimizeCbModel(mtmp,'max','one'); 
    WTFlux(i,:) = (sol1n.x)';
  end
  potentialFail = abs(WTFlux(:,[536,1577])) < 1e-7;
  pFsum = sum(potentialFail(:,1) & potentialFail(:,2));
  if pFsum > 0
    error(strcat(num2str(pFsum)),' unlikely solutions encountered from geoFBA')
  end
end



