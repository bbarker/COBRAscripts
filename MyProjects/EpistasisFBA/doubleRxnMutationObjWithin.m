function [grMat, epiMult, smfAdj] = doubleRxnMutationObjWithin(modelOld, model, atime, mutMat, WTflux, SPvec, lMOMA)

n = length(mutMat);
grMat = zeros(n,n);
epiMult = zeros(n,n);
smfAdj = zeros(1,n);

%fluxList = zeros(2^n,length(model.lb));
nrxns = length(model.rxns);

solOld = optimizeCbModel(modelOld,'max',1); %use L2 norm
solFBA = optimizeCbModel(model,'max');

disp('Check fitnesses old:');
disp(solOld.f);

if nargin < 5
  WTflux = solOld.x;
end
if nargin < 6
  SPvec = ones(1,length(modelWT.lb));
end
if nargin < 7
  lMOMA = true;
end

if lMOMA
  %Change f and g!
  [solAda, solutionWT] = linearMOMAOneShot(model,WTflux,atime,0.02,.07);
else
  [solAda,solutionWT,totalFluxDiff,solStatus] = MOMATest(modelOld,model,'max',false,true,atime,WTflux,SPvec);
end

disp('Check fitnesses2:');
disp([solOld.f solutionWT.f solAda.f]);

grRow = zeros(1,n);

parfor i = 1:n
  grRow = zeros(1,n);
  for j = i+1:n
    if mutMat(i,3) ~= mutMat(j,3)
      wvec = zeros(1,nrxns);
      wvec(mutMat(i,3)) = mutMat(i,2);
      wvec(mutMat(j,3)) = mutMat(j,2);
      %Try this to reduce double mutant fitness going to 0:
      %wvec = wvec/2; %this works well, now try geometric
      %wvec(mutMat(i,3)) = sign(mutMat(i,2))*sqrt(abs(mutMat(i,2)));
      %wvec(mutMat(j,3)) = sign(mutMat(j,2))*sqrt(abs(mutMat(j,2)));
      if lMOMA
	[solKO, solutionWT] = linearMOMAOneShot(model,WTflux,atime,0.02,.07,wvec);
      else
	[solKO,solutionWT,totalFluxDiff,solStatus] = MOMATest(modelOld,model,'max',false,true,atime,WTflux,SPvec,wvec); %WTflux or solAda.x?
      end
      grRow(j) = solKO.f;
    end
  end
  grMat(i,:) = grRow; 
end

for i = 1:n
  wvec = zeros(1,nrxns);
  wvec(mutMat(i,3)) = mutMat(i,2);
  %wvec(mutMat(i,3)) = mutMat(i,2)/2;
  if lMOMA
    [solKO, solutionWT] = linearMOMAOneShot(model,WTflux,atime,0.02,.07,wvec);
  else
    [solKO,solutionWT,totalFluxDiff,solStatus] = MOMATest(modelOld,model,'max',false,true,atime,WTflux,SPvec,wvec); %WTflux or solAda.x?
  end
  smfAdj(i) = solKO.f;
end

for i=1:n
  for j=i+1:n
    if mutMat(i,3) ~= mutMat(j,3)
      epiMult(i,j) = grMat(i,j)/solAda.f - smfAdj(i)*smfAdj(j)/(solAda.f^2); 
    else
      epiMult(i,j) = 0;
    end
  end
end

