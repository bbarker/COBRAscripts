function [grTens, rxnTens] = randomBeneficialPsets(simprefix, modelOld, model, atime, SPvec, WTflux, mutMat, rsetsize, nsamp, beneThresh, lMOMA)

if nargin < 10
  beneThresh = 0; %Threshold for relative fitness to sample mutations from.
end
numFound = sum(mutMat(:,:,1)'>=(beneThresh));


if nargin < 11
  lMOMA = true;
end

geneRxns = any(model.rxnGeneMat,2);
geneRxns = find(geneRxns);
geneRxns = geneRxns(find(numFound));
nGRxns = length(geneRxns);
mutMat = squeeze(mutMat(find(numFound),:,:));
numFound = numFound(find(numFound));

mutSZ = size(mutMat);
maxSMuts = mutSZ(2);

%grRateListRand10_3 = adaptiveMomaPSetMutation('grRateListRand10_3', mm, mmEth, mutMat3(rand10_3,3), 1680, SPvec, mutMat3(rand10_3,2),WTfluxEthMod);

%[grRateList, rxnlistrev] = adaptiveMomaPSetMutation(simprefix, modelOld, model, rxnlist, atime, SPvec, MWeight, WTflux)
%rand10_4=(1043-57)+randNoResample(57,10)
msz = size(mutMat);
% mutMat is (fitness, flux restriction, reaction #)
grTens = zeros(nsamp,2^rsetsize);
rxnTens = zeros(nsamp,rsetsize);

zvec = zeros(1,length(model.lb));
if lMOMA
  %Change f and g!
  [solAda, solutionWT] = linearMOMAOneShot(model,WTflux,atime,0,1,zvec);
else %fix later if needed:
  [solAda,solutionWT,totalFluxDiff,solStatus] = MOMATest(modelOld,model,'max',false,true,atime,WTflux,SPvec);
end

%nbene = sum(mutMat(:,1)/solAda.f > beneThresh);
%This needs to be altered to deal with new mutMat structure.
%Assume all are beneficial, or at least of interest.
failednum = 0;
successes = 0;
for i = 1:nsamp
  randrxns = zeros(rsetsize,3);
  succ = 0;
  randGRxns = [];
  randFluxIdx = zeros(rsetsize,1);
  while (succ == 0)
    %Probably need to parallelize the testing as well
    %Repeatedly sampling fluxes from the same random reaction
    %set may cause problems if there are troublesom reactions.
    modelMut = model;
    randGRxns = randNoResample(nGRxns, rsetsize);
    for j = 1:rsetsize
      randFluxIdx(j) = randi(min(numFound(randGRxns(j)),maxSMuts),1);
      modelMut.lb(geneRxns(randGRxns(j))) = mutMat(randGRxns(j),randFluxIdx(j),2);
      modelMut.ub(geneRxns(randGRxns(j))) = modelMut.lb(geneRxns(randGRxns(j)));
    end
    if lMOMA
      [solKO, solutionWT] = linearMOMAOneShot(modelMut,WTflux,atime,0,1,zvec);
    else
      [solKO,solutionWT,totalFluxDiff,solStatus] = MOMATest(modelOld,modelMut,'max',false,true,atime,WTflux,SPvec,wvec); %WTflux or solAda.x?
    end
    if solKO.stat > 0.5
      succ = 1;
      successes = successes + 1;
    else
      failednum = failednum + 1;
    end
  end
  for j = 1:rsetsize
    rGR = randGRxns(j);
    rxnNum = geneRxns(rGR);
    %rsmut = randi(min(numFound(rGR),maxSMuts),1);
    rsmut = randFluxIdx(j);
    rflux = mutMat(rGR,rsmut,2);
    rfit = mutMat(rGR,rsmut,1);
    randrxns(j,:) = [rfit rflux rxnNum];
  end  
  [grRateList, rxnlistrev] = beneficialMomaPSetMutation('grRateList', modelOld, model, randrxns, atime, SPvec ,WTflux);
  grTens(i,:) = grRateList;
  rxnTens(i,:) = rxnlistrev;
end
dlmwrite(strcat(simprefix,'_grRateMat.csv'), grTens, ',');
dlmwrite(strcat(simprefix,'_rxnlist.csv'), rxnTens, ',');
disp(failednum);
disp(successes);