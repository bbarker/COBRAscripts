function [grTens, rxnTens] = randomAdaptivePsets(simprefix, modelOld, model, atime, SPvec, WTflux, mutMat, rsetsize, nsamp, beneThresh, lMOMA)

if nargin < 10
  beneThresh = 1.05; %Threshold for relative fitness to sample mutations from.
end
if nargin < 11
  lMOMA = true;
end
%grRateListRand10_3 = adaptiveMomaPSetMutation('grRateListRand10_3', mm, mmEth, mutMat3(rand10_3,3), 1680, SPvec, mutMat3(rand10_3,2),WTfluxEthMod);

%[grRateList, rxnlistrev] = adaptiveMomaPSetMutation(simprefix, modelOld, model, rxnlist, atime, SPvec, MWeight, WTflux)
%rand10_4=(1043-57)+randNoResample(57,10)
msz = size(mutMat);

grTens = zeros(nsamp,2^rsetsize);
rxnTens = zeros(nsamp,rsetsize);
if lMOMA
  %Change f and g!
  [solAda, solutionWT] = linearMOMAOneShot(model,WTflux,atime,0.02,.07);
else
  [solAda,solutionWT,totalFluxDiff,solStatus] = MOMATest(modelOld,model,'max',false,true,atime,WTflux,SPvec);
end
nbene = sum(mutMat(:,1)/solAda.f >= beneThresh);

for i = 1:nsamp
  randrxns = (msz(1)-nbene) + randNoResample(nbene, rsetsize);
  [grRateList, rxnlistrev] = adaptiveMomaPSetMutation('grRateList', modelOld, model, mutMat(randrxns,:), atime, SPvec ,WTflux);
  grTens(i,:) = grRateList;
  rxnTens(i,:) = rxnlistrev;
end
dlmwrite(strcat(simprefix,'_grRateMat.csv'), grTens, ',');
dlmwrite(strcat(simprefix,'_rxnlist.csv'), rxnTens, ',');
