function epiDat = randomEpistasisSampler(model, InitFlux, smutDat, numFound, beneOnly)
sthresh = 0.001;
maxEpi = 4000000
growthRxn = find(model.c);
zthresh = 1e-8;
objWeight = 0;



epiDat = zeros(maxEpi,9);
%Fields are (rxn1,rxn2,flux1,flux2,rfit1,rfit2,rfitDbl,epsilon,dblStat);

%EpiInfo = struct('epiCount',0,'intTests',0,'posEpi',0,'NegEpi',0);

if beneOnly
  numFound = sum(smutDat(:,:,1)'>=(1+sthresh));
end

nRxns = length(InitFlux);
geneRxns = any(model.rxnGeneMat,2);
geneRxns = find(geneRxns);
geneRxns = geneRxns(find(numFound));
nGRxns = length(geneRxns);
smutDat = squeeze(smutDat(find(numFound),:,:));
numFound = numFound(find(numFound));

smutSZ = size(smutDat);
maxSMuts = smutSZ(2);

zvec = zeros(1,nRxns);


%Sampling method: equal chance of picking a reaction (biased
%towards reactions with fewer beneficial mutations).
parfor i = 1:maxEpi
  rGR1 = randi(nGRxns,1);
  rGR2 = randNoResample(nGRxns, 1, rGR1);
  rxnNum1 = geneRxns(rGR1);
  rxnNum2 = geneRxns(rGR2);
  rsmut1 = randi(min(numFound(rGR1),maxSMuts),1);
  rsmut2 = randi(min(numFound(rGR2),maxSMuts),1);
  rflux1 = smutDat(rGR1,rsmut1,2);
  rflux2 = smutDat(rGR2,rsmut2,2);
  rfit1 = smutDat(rGR1,rsmut1,1);
  rfit2 = smutDat(rGR2,rsmut2,1);
  modelMut = model;
  modelMut.lb(rxnNum1) = rflux1
  modelMut.ub(rxnNum1) = rflux1;
  modelMut.lb(rxnNum2) = rflux2
  modelMut.ub(rxnNum2) = rflux2;

  [solutionDel,solutionWT] = linearMOMAOneShot(modelMut,InitFlux,objWeight,0,1,zvec);
  rfitDbl = solutionDel.f/InitFlux(growthRxn);
  epi = rfitDbl - rfit1*rfit2;
  if epi > zthresh
    disp([epi rxnNum1 rxnNum2]);
  end
  epiDat(i,:) = [rxnNum1 rxnNum2 rflux1 rflux2 rfit1 rfit2 rfitDbl epi solutionDel.stat];
end