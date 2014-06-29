function modelOut = makeKishonyDeLunaModel(modelIn)

%Growth rate should be similar to standard YPD medium.

modelOut = modelIn;


% Some constants
LBmin = -1000;


%
% Correct to appropriate genetic background (TODO!)
%

% Adjust glucose uptake:
idx = find(strcmp(modelOut.rxnNames, 'D-glucose exchange'));
modelOut.lb(idx) = -10;

% Change nitrogen source:
idx = find(strcmp(modelOut.rxnNames, 'ammonium exchange'));
modelOut.lb(idx) = 0;
idx = find(strcmp(modelOut.rxnNames, 'L-proline exchange'));
modelOut.lb(idx) = -1;

%
% Verify uptakes
%
[selExc, selUpt] = findExcRxns(modelOut);
modelOut.rxnNames(selUpt)
