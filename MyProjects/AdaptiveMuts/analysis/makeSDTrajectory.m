function [AdaptedFitness, RelativeFitness] = makeSDTrajectory(modelOld,modelNew,WTflux,initW,finalW,step,simname,SPvec, genenums)

%Maybe make a parallel version that calculates mutants at each step
%For plotting a figure in a paper.
FluxThresh = .05;
minCut = 1e-4;

entries = [initW:step:finalW];
if nargin < 8
  SPvec=ShadowPricesHeuristic(modelOld,WTflux);
end
%The last argument of the line below may need to be tweaked.
%minL2New = optimizeCbModel(modelNew,'max',0.1);
nentries = length(entries);
%optOld = optimizeCbModel(modelOld);
%optNew = optimizeCbModel(modelNew);
%solAdaPrev.x = optNew.f*WTflux/WTflux(find(modelOld.c)); %This is dubious, so ignore first entry for now.
glist = {''};
if nargin < 9
  glist = [glist; modelNew.genes];
else
  glist = [glist; modelNew.genes(genenums)];
end
AdaptedFitness = zeros(length(glist), length([initW:step:finalW]))'; 
for j=1:length(glist)
  [modelDel,hasEffect,constrRxnNames] = deleteModelGenes(modelNew,glist{j}); %uses modified script
  parfor i = 1:nentries;
    [solAda,solutionWT,totalFluxDiff,solStatus] = MOMATest(modelOld,modelDel,'max',false,true,entries(i),WTflux,SPvec);
    AdaptedFitness(i,j) = solAda.f;
    disp(simname);
    disp([entries(i) solAda.f solutionWT.f]);
    %Why this does not work in parallel sometimes?
  end
end

RelativeFitness = zeros(length(glist)-1, length([initW:step:finalW]))';
for j = 1:(length(glist)-1)
  for i = 1:length([initW:step:finalW])
    RelativeFitness(i,j) = AdaptedFitness(i,j+1) / AdaptedFitness(i,1);
  end
end
save(strcat(simname,'.mat'), 'initW', 'finalW', 'step', 'AdaptedFitness', 'RelativeFitness');

