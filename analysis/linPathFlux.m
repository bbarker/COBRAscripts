function [qVecStart, qVecPrior] linPathFlux(m, varargin)
% Measures attribution of a (metabolite) flux along a linear
% pathway in terms of the start of the pathway. 
%
% metNames are used to define pathways as these are the most
% general possibility and can even encompass multiple compartments; 
% for a more specific use, reactions can be used (see isRxn).

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
% COBRA model
%
p.addRequired('linPath', @IPcheck_linPath);
% A cell array of metNames or alternatively may contain
% some rxn ids (see isRxn).
% 
%
%%% Optional INPUT %%%
%
p.addParamValue('isRxn', [false], @IPcheck_boolVec);
% A Boolean vector of the same length as linPath
% with true entries corresponding to linPath entries
% that are rxn ids instead of metNames 
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


end % of linPathFlux


%%%%%%%%%%%%%     Input Parser Functions     %%%%%%%%%%%%%

function TF = IPcheck_linPath(x)
TF = true;
%check x is a cell array
if ~iscell(x)
    error('linPath must be cell array of strings or cells of strings');
end
% check x is 1-d
sz_x = size(x);
if min(sz_x) ~= 1 || numel(sz_x) > 2
    error('x is not one dimensional');
end
%check each entry is a string or a cell
%actually, save this as a model-based check later
end % of IPcheck_linPath

function TF = IPcheck_boolVec(x)
TF = true;
if ~islogical(x) || ~isvector(x)
    error('input is not a logical vector');
end
end % of IPcheck_boolVec

%%%%%%%%%%%%%   End Input Parser Functions   %%%%%%%%%%%%%
