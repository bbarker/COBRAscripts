function [transRxns,nonTransRxns] = findTransRxns2(model,inclExc, ... 
    rxnInds,inclObjAsExc,irrevFlag)
% like findTransRxns in the COBRAtoolbox, but doesn't rely on 
% compartment labels. Instead, it relies on the metNames field being
% a unique species identifier.
%
%findTransRxns identify all transport reactions in a model, which are
%defined as reactions involved with metabolites in more than 1 compartment
%
%INPUT
% model             COBRA model structure
%
%OPTIONAL INPUT
% inclExc           include exchange reactions as transport?
%                   (Default = false)
% rxnInds           indices of reactions to test for transport activity.
%                   (default = test all columns of model.S)
% inclObjAsExc      include objective as an exchange reaction? this is
%                   passed to findExcRxns. (default = false)
% irrevFlag         is model irreversible? this is passed to findExcRxns.
%                   (default=false)
%
%OUTPUT
% transRxns         transport reactions in rxnInds
% nonTransRxns      non-transport reactions in rxnInds

if nargin < 2
    inclExc = false;
end
if nargin < 3
    rxnInds = 1:size(model.S, 2);
end
if nargin < 4
    inclObjAsExc = false;
end
if nargin < 5
    irrevFlag = false;
end

if inclExc
    % findExcRxns returns boolean vector
    isExc0 = findExcRxns(model,inclObjAsExc,irrevFlag);
    % subset to rxnInds
    isExc = isExc0(rxnInds);
else
    isExc = zeros(length(rxnInds), 1);
end

%boundaryRxns = find(  sum(boolean(Sabs)) == 1  );
%istrans(boundaryRxns) = 1;


function istrans = isTransport(model, rxnNum, metName)
Sabs = abs(model.S);
nmrxns = length(rxnNum);
istrans = zeros(1,nmrxns);
for i = 1:nmrxns
    rxnmets = find(model.S(:, rxnNum(i)));
    rxnMetNames = model.metNames(rxnmets);
    nameMatches = find(strcmp(rxnMetNames, metName));
    if length(nameMatches) > 1
        istrans(i) = 1;
    end
end
end % End of isTransport
