function [residSum, min2normDist, fitnessVec, WeightsGenes, DecGenes, IncGenes, RxnDeltaList] = makeTrajectory(modelOld,modelNew,WTflux,initW,finalW,step,simname,SPvec,lMOMA,fw,gw)

%Maybe make a parallel version that calculates mutants at each step
%For plotting a figure in a paper.
FluxThresh = .01;
minCut = 1e-4;

entries = [initW:step:finalW];
if nargin < 8
  %SPvec=ShadowPricesHeuristic(modelOld,WTflux);
  SPvec = ones(1,length(modelNew.lb));
end
if nargin < 9
  lMOMA = true;
end
if nargin < 10
  fw = 0.02;
end
if nargin < 11
  gw = 0.07
end

minL2New = optimizeCbModel(modelNew,'max',2);
nentries = length(entries);
%optOld = optimizeCbModel(modelOld);
optNew = optimizeCbModel(modelNew);
residSum = zeros(1,nentries);
min2normDist = zeros(1,nentries);
fitnessVec = zeros(1,nentries);
WeightsGenes = [];
DecGenes = [];
IncGenes = [];
RxnDeltaList = [];
RxnDiffPrev = modelNew.c*0;
RxnDelta = modelNew.c*0;
solAdaPrev.x = optNew.f*WTflux/WTflux(find(modelOld.c)); %This is dubious, so ignore first entry for now.
solAdaPrev.f=-1;
for i = 1:nentries;
  if lMOMA
    %Change f and g!
    [solAda, solutionWT] = linearMOMAOneShot(modelNew,WTflux,entries(i),fw,gw);
  else
    [solAda,solutionWT,totalFluxDiff,solStatus] = MOMATest(modelOld,modelNew,'max',false,true,entries(i),WTflux,SPvec);
  end
  RxnDiff = abs(solAda.x-WTflux) > max(minCut,abs(FluxThresh*WTflux));
  if sum(RxnDiff ~= RxnDiffPrev) > 0
    WeightsGenes = [WeightsGenes [entries(i); modelNew.rxnGeneMat' * RxnDiff]];
  end
  RxnDelta = abs(solAda.x/solAda.f) - abs(solAdaPrev.x/solAdaPrev.f);
  IncGenes = [IncGenes modelNew.rxnGeneMat' * (RxnDelta > max(minCut,abs(FluxThresh*solAdaPrev.x/solAdaPrev.f)))];
  DecGenes = [DecGenes modelNew.rxnGeneMat' * (RxnDelta < -max(minCut,abs(FluxThresh*solAdaPrev.x/solAdaPrev.f)))];
  RxnDeltaList = [RxnDeltaList RxnDelta];
  residSum(i) = norm(solAda.x - WTflux);
  min2normDist(i) = norm(solAda.x) - norm(minL2New.x); 
  fitnessVec(i) = solAda.f;
  RxnDiffPrev = RxnDiff;
  solAdaPrev = solAda; 
end
save(strcat(simname,'.mat'), 'residSum', 'min2normDist', 'fitnessVec', 'WeightsGenes', 'RxnDeltaList', 'DecGenes', 'IncGenes');
