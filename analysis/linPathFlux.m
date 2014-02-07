function [qVecStart, pathInfo] = linPathFlux(m, varargin)
% TODO:
% Remove transport and exchange rxns from analysis
%
%
% Measures attribution of a (metabolite) flux along a linear
% pathway in terms of the start of the pathway. 
%
% There are two ways to look at this:
% 1) Percent of X converted to Y along pathway X...Y
% 2) Percent of Y derived from X along pathway X...Y
% For now we focus only on 1.
%
% metNames are used to define pathways as these are the most
% general possibility and can even encompass multiple compartments; 
% for a more specific use, reactions can be used..

p = inputParser;
p.FunctionName = 'linPathFlux';
p.StructExpand = false;
p.CaseSensitive = true;
if isfield(p, 'PartialMatching')
    p.PartialMatching = false;
end

%%% Required INPUT %%%
%
p.addRequired('m');
% COBRA *irreversible* model
%
p.addRequired('linPath', @IPcheck_linPath);
% A cell array of metNames or alternatively may contain
% some rxn ids as a cell array of rxn ids, where
% 
p.addRequired('flux', @IPcheck_NumVec);
% irreversible model flux vector
%
%%% Optional INPUT %%%
%
%
p.addParamValue('rxns', {}, @IPcheck_rxns);
%
% OUTPUT
%
%qVecStart
% The ratio of flux at each point in the pathway that 
% can be attributed to the initial (metabolite) flux
%
%qVecPrior
% The ratio of flux at each point in the pathway that 
% can be attributed to the previous flux in the pathway


p.parse(m, varargin{:});
arg = p.Results;

linPathRxns = cell(1, length(arg.linPath)-1);
linPathSubs = cell(1, length(arg.linPath)-1);
allOutRxns  = cell(1, length(arg.linPath)-1);
allOutSubs  = cell(1, length(arg.linPath)-1);
qVecStart = zeros(1, length(arg.linPath)-1);
qVecPrior = zeros(1, length(arg.linPath)-1);
pathLen = length(arg.linPath);
%Convert linPath into reaction sets.
%Also list alternative branch points and |substrate coeffs|.
for i = 1:(pathLen-1)
    %Get all substrates
    substrates = find(strcmp(arg.m.metNames, arg.linPath{i}));
    products   = find(strcmp(arg.m.metNames, arg.linPath{i+1}));
    %Get all forward reactions
    for j = 1:length(substrates)
        s = substrates(j);
        AsRxns = find((arg.m.S(s, :) < 0));
        allOutRxns{i} = [allOutRxns{i} AsRxns];
        allOutSubs{i} = [allOutSubs{i} full(-arg.m.S(s, AsRxns))];
        spRxns = [];
        for k = 1:length(products)
            p = products(k);
            spRxns = [spRxns find((arg.m.S(s, :) < 0) & (arg.m.S(p, :) > 0))];
        end
        linPathRxns{i} = [linPathRxns{i} spRxns];
        linPathSubs{i} = [linPathSubs{i} full(-arg.m.S(s, spRxns))];
    end
end

%Optionally, if reaction sets are explicitly listed for any
%pathway components, override with them here:
if length(arg.rxns) > 0
    for i = 1:(pathLen-1)
	linPathRxns{i} = zeros(size(arg.rxns{i}));
	for j = 1:length(linPathRxns{i})
            if length(arg.rxns{i}(j)) > 0
	        linPathRxns{i}(j) = ...
                    find(strcmp(arg.m.metNames, arg.rxns{i}(j)));
            end
	end
    end
end

% For each step in the pathway, compute the 
% proportion of metabolite flux going in to the 
% next metabolite.

for i = 1:(pathLen-1)    
    % We need to account for substrate stoichiometry coefficients
    % when we are looking for amount of X converted to Y.
    idx = linPathRxns{i}
    adx = allOutRxns{i} 
    is = linPathSubs{i}
    as = allOutSubs{i}
    qVecPrior(i) = (arg.flux(idx)' * is') / (arg.flux(adx)' * as')
    if i > 1
        qVecStart(i) = qVecPrior(i) * qVecStart(i-1); 
    else
        qVecStart(i) = qVecPrior(i);
    end 
end

pathInfo.linPathRxns = linPathRxns;
pathInfo.linPathSubs = linPathSubs;
pathInfo.allOutRxns = allOutRxns;
pathInfo.allOutSubs = allOutSubs;
pathInfo.qVecStart = qVecStart;
pathInfo.qVecPrior = qVecPrior;
end % of linPathFlux

%%%%%%%%%%%%%     Input Parser Functions     %%%%%%%%%%%%%

function TF = IPcheck_linPath(x)
TF = true;
%check x is a cell array
if ~iscell(x)
    error('linPath must be a cell array');
end
% check x is 1-d
sz_x = size(x);
if min(sz_x) ~= 1 || numel(sz_x) > 2
    error('x is not one dimensional');
end
%check each entry is a string
for i = 1:length(x)
    if ~isstr(x{i})
        error(['entry ' num2str(i) ' not a string.']); 
    end
end
end % of IPcheck_linPath

function TF = IPcheck_boolVec(x)
TF = true;
if ~islogical(x) || ~isvector(x)
    error('input is not a logical vector');
end
end % of IPcheck_boolVec

function TF = IPcheck_NumVec(x)
TF = true;
if ~isvector(x)
    error('flux must be a vector.')
elseif ~isnumeric(x)
    disp(x)
    error('flux must be numeric.')
end
end % end of IPcheck_NumVec

%%%%%%%%%%%%%   End Input Parser Functions   %%%%%%%%%%%%%
