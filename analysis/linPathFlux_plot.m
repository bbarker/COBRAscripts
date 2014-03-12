function linPathFlux_plot(qMat, baseColor, xVec, xLabel, pathIdxs, linPath)
% Given a consecutive list of outputps from a script calling linPathFlux 
% (in matrix form: qMat), plots them in the order presented as 
% line plots, with the first entries the darkest shade of the specified
% color
%
%INPUT
%
% qMat:   array of outputs from linPathFlux

%Optional INPUTS
%
% linPath: The same path passed to linPathFlux.
%          Used to construct legends.
%
% pathIdxs: Only plot these reactions.
%

white = [1 1 1];

[nIter, pathLenM1] = size(qMat);

if ~exist('baseColor', 'var')
    baseColor = [0 0 0]; % black;
end
maxBC = max(baseColor);

if ~exist('xVec', 'var')
    xVec = 1:nIter;
end 

if ~exist('xLabel', 'var')
    xLabel = '';
end 

if ~exist('pathIdxs', 'var')
    pathIdxs = 1:pathLenM1;
end 

figure();
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
set(gca, 'FontSize', 28);
hold all;

for i = pathLenM1:-1:1
    if any(ismember(pathIdxs, i))
        iColor = baseColor + ((i-1)/pathLenM1) * (1-maxBC) * white;
        plot(xVec, qMat(:, i)', 'Color', iColor, 'LineWidth', 3);
    end
end
xlabel(xLabel);
if exist('linPath', 'var')
    ylabel(['Fraction of ' linPath{1} ' converted.']);
else
    ylabel('Fraction of entry metabolite converted.');
end

if exist('linPath', 'var')
    legEntries = cell(1, pathLenM1);
    for i = 1:pathLenM1
        if any(ismember(pathIdxs, i))
            legEntries{i} = [linPath{i} ' -> ' linPath{i+1}]; 
        end
    end
    legend(legEntries{:});
end

