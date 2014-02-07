function [qVecStart, qVecPrior] = linPathFlux(m, varargin)
qVecStart = 0;
qVecPrior = 0;
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
arg = p.Results

%Convert linPath into reaction sets.
linPathRxns = cell(1, length(arg.linPath)-1);
pathLen = length(arg.linPath);
for i = 1:(pathLen-1)
    %Get all substrates
    substrates = find(strcmp(arg.m.metNames, arg.linPath{i}));
    products   = find(strcmp(arg.m.metNames, arg.linPath{i+1}));
    %Get all forward reactions
    for j = 1:length(substrates)
        s = substrates(j);
        for k = 1:length(products)
            p = products(k);
            linPathRxns{i} = union(linPathRxns{i}, ...
                find((arg.m.S(s, :) < 0) & (arg.m.S(p, :) > 0)));
        end
    end
    linPathRxns{i} = setdiff(linPathRxns{i}, [0]);
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

%%%%%%%%%%%%%   End Input Parser Functions   %%%%%%%%%%%%%
