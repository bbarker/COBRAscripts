function [grRateList, rxnlistrev] = beneficialMomaPSetMutation(simprefix, modelOld, model, mutMat, atime, SPvec, WTflux, lMOMA)
%How do we deal with crashes? Answer: extend this later to
%take its own output and extend it with additional genes. For
%large n, may need to modify the program to add a superficial
%loop, instead of just using a parfor, thus allowing 
%incremental saving.

zvec = zeros(1,length(model.lb));

rxnlist = mutMat(:,3);
MWeight = mutMat(:,2); %actually restriction here
rxnlistrev = fliplr(rxnlist);
n = length(rxnlistrev);
psmax = -1+2^n;

grRateList = zeros(1,2^n);
%fluxList = zeros(2^n,length(model.lb));
nrxns = length(model.rxns);

solOld = optimizeCbModel(modelOld,'max',1); %use L2 norm
solFBA = optimizeCbModel(model,'max');

if nargin < 7
  WTflux = solOld.x;
  WTflux(find(model.c)) = solFBA.x(find(model.c));
end
if nargin < 8
  lMOMA=true;
end

if lMOMA
  %Change f and g!
  [solAda, solutionWT] = linearMOMAOneShot(model,WTflux,atime,0,1,zvec);
else %fix later if needed:
  [solAda,solutionWT,totalFluxDiff,solStatus] = MOMATest(modelOld,model,'max',false,true,atime,WTflux,SPvec);
end

disp('Check fitnesses for WT are the same:');
disp([WTflux(find(model.c)) solAda.x(find(model.c))]);

disp('Check fitnesses2:');
disp([solOld.f solutionWT.f solAda.f]);


parfor i = 0:psmax
  modelMut = model;
  bstrvec = boolean(str2numvec(dec2bin(i,n)));
  subrxnlist = rxnlist(bstrvec);
  %wvec = zeros(1,nrxns);
  svec = zeros(1,length(subrxnlist));
  for j = 1:length(subrxnlist)
    listidx = find(ismember(rxnlist,subrxnlist(j)));
    %if length(listidx) > 1
    %  disp(listidx);
    %  disp(j);
    %  disp(rxnlist);
    %  disp(subrxnlist);
    %end
    disp('ListIdx dbg:');
    disp(listidx);
    disp(rxnlist(listidx));
    svec(j) = mutMat(listidx,1)/solAda.f - 1;
  end
  for j = 1:length(subrxnlist)
    listidx = find(ismember(rxnlist,subrxnlist(j)));
    %wvec(subrxnlist(j)) = MWeight(listidx);
    modelMut.lb(subrxnlist(j)) = MWeight(listidx);
    modelMut.ub(subrxnlist(j)) = MWeight(listidx);
    %Also consider trying sum instead of max, or some variation of these.
  end
  %if length(subrxnlist) == 1
  %  disp('Before and after');
  %  disp(wvec(subrxnlist(j)));
  %end
  %if length(subrxnlist) == 1
  %  disp(wvec(find(wvec)));
  %end
  if lMOMA
    [solKO, solutionWT] = linearMOMAOneShot(modelMut,WTflux,atime,0,1,zvec);
  else
    [solKO,solutionWT,totalFluxDiff,solStatus] = MOMATest(modelOld,modelMut,'max',false,true,atime,WTflux,SPvec,wvec); %WTflux or solAda.x?
  end
  grRateList(i+1) = solKO.f;
  if length(subrxnlist) == 1
    listidx = find(ismember(rxnlist,subrxnlist));
    disp([mutMat(listidx,:) solKO.f] );
  end
  
%  fluxList(i+1,:) = solKO.x;
end



