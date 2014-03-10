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


%I f model includes compartment identifiers, remove these:
for i = 1:length(model.metNames)
    %Yeast 7 style:
    model.metNames{i} = regexprep(model.metNames{i}, ...
        '(\s+)?\[.+](\s+)?$', '');
end

%boundaryRxns = find(  sum(boolean(Sabs)) == 1  );
% use isExc instead for now
%
isNonexchTrans = isTransport(model, rxnInds);


% get rxn abbreviations for all rxns in rxnInds
rxnAbbrevs=model.rxns(rxnInds);
% if inclExc==1, exchange rxns will have isExc==1, and should be counted as
% transport rxns; else, all isExc will be 0.
transRxns = rxnAbbrevs(isNonexchTrans==1 | isExc==1);
nonTransRxns = setdiff(rxnAbbrevs, transRxns);

end % end of findTransRxns2

function istrans = isTransport(model, rxnNum)
Sabs = abs(model.S);
nmrxns = length(rxnNum);
istrans = zeros(nmrxns, 1);
for i = 1:nmrxns
    rxnmets = find(model.S(:, rxnNum(i)));
    rxnMetNames = model.metNames(rxnmets);
    for j = 1:length(rxnMetNames)
        metName = rxnMetNames{j};
        nameMatches = find(strcmp(rxnMetNames, metName));
        if length(nameMatches) > 1
            istrans(i) = 1;
            continue;
        end
    end
end
end % End of isTransport
