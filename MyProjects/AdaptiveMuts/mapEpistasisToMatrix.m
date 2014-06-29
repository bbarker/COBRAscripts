function [hmMat, posPerc, hmCnt, fminmax] = mapEpistasisToMatrix(epiDat, matlen)
zthresh = 1e-8;
eCut = 0.01;
hmMat = zeros(matlen, matlen);
hmCnt = zeros(matlen, matlen);
posPerc = zeros(matlen, matlen);

%epiDat fields are (rxn1,rxn2,flux1,flux2,rfit1,rfit2,rfitDbl,epsilon);

minSMF = min(flat(epiDat(:,[5 6]))) 
maxSMF = max(flat(epiDat(:,[5 6])))
SMF_inc = (maxSMF - minSMF)/matlen
fminmax = [minSMF maxSMF];

epiSZ = size(epiDat);
errterm = 1e-7;
for i = 1:epiSZ(1)
  x = ceil((epiDat(i,5)-minSMF+errterm)/SMF_inc);
  y = ceil((epiDat(i,6)-minSMF+errterm)/SMF_inc);
  if y == matlen+1
    y = matlen;
  end
  if x == matlen+1
    x = matlen;
  end
  hmMat(x,y) = hmMat(x,y) + epiDat(i,8);
  hmCnt(x,y) = hmCnt(x,y) + 1;
  if epiDat(i,8) > eCut
    posPerc(x,y) = posPerc(x,y) + 1;
  end
end

hmMat = hmMat ./ hmCnt;
posPerc = 100*posPerc ./ hmCnt; 
%naninf = isnan(A) + isinf(A);
naninf = isnan(hmMat); %only nans should be an issue
hmMat(boolean(naninf)) = 0;
naninf = isnan(posPerc); %only nans should be an issue
posPerc(boolean(naninf)) = 0;
