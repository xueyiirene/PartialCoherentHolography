function [pts, permidx] = function_poissonDisc_fixgrid(sGrid,spacing)

% Purpose:
% N-dimensional poisson disc sampling function. This can also be used to
% randomly sample k pts from N-dimensional space with a minimum separation
% distance.
%
% Inputs:
% sizeI -   [required] Size of volume from which points are to be 
%           sampled
% spacing - [required] Minimum sepration distance between points
% nPts -    [Default is 0] if nPts = 0 For poisson disc sampling.
%           nPts = k to sample k-pts from N-dimensional space with
%           minimum separation distance
% showIter - [Default is 0] If showIter == 1, this option can be used to 
%            see how points are generated through each iteration. It can be
%            useful when code is taking a long time to generate points. 
%
% Output:
% pts - All eligible points
%
%
% Example:
% 1. Poisson disc sampling in 2-dimensional space. 
% sizeI = [512,512];
% spacing = 30;
% pts = poissonDisc(sizeI,spacing);
%
% 2. Sample k-pts in 3-dimensional space.
% sizeI = [512,512,192];
% spacing = 6;
% nPts = 10000;
% pts = poissonDisc(sizeI,spacing,nPts);
% 
% 3. Show iteration progress from poisson disc sampling in 2-dimension
% sizeI = [512,512];
% spacing = 6;
% nPts = 0;
% showIter = 1;
% pts = poissonDisc(sizeI,spacing,nPts,showIter);

% Mohak Patel, Brown University, 2016

%%%%%%% Initial parameters setup
% Parsing inputs and setting default values
if nargin == 3;  end
if nargin == 2;  nPts = 0; end

% Setting properties for iterations
ndim = 2;   % Number of Dimensions
k = 5;  % Number of 'dart' tries in each grid.
dartFactor = 4; %Select number of sample data in each iterations. Change it to
% reduce run time for code. Have to play around with number. 


%%%%%%% Making Grid read for iterations
%Make grid size such that there is just one pt in each grid
dm = spacing/sqrt(ndim);    % grize cell size [Bridson 2007]

% Arrays to show eligible grids for dart throws and keeping score of darts
% thrown in a particular grid
emptyGrid = logical(ones(size(sGrid,1),1)); %Eligible Grids
nEmptyGrid = sum(emptyGrid);    %Number of eligible Grids
scoreGrid = zeros(size(emptyGrid)); %Score of darts thrown in Grid

% Darts to be thrown per iterations
% This hugely influences speed of the algorithm. Change dartFactor for it. 
nPts=size(sGrid,1);
ndarts = round(nPts/dartFactor);

%%%%%%%%% Iterative process to generate points
% Initialize parameters
ptsCreated = 0;
pts = [];
iter = 0;

% Start Iterative process
while ptsCreated<nPts && nEmptyGrid >0   
    %Thrown darts in eligible grids
    availGrid = find(emptyGrid == 1);   %Eligible grids for dart throw
    dataPts = min([nEmptyGrid,ndarts]); % Darts to be thrown
    p = datasample(availGrid,dataPts,'Replace',false); %Select grids for darts
    tempPts = sGrid(p,:); %Dart throw!!!
    
    % Find good dart throws
    [~,D] = knnsearch([pts;tempPts],tempPts,'k',2); %Finding distance between all darts(pts)
    D = D(:,2); 

    eligiblePts = D>spacing; %elgible pts should also have minimum separation distance
    
    scorePts = tempPts(~eligiblePts,:); %Keep score from bad dart throws :(
    tempPts = tempPts(eligiblePts,:);   % Save good dart throws :)
    if sum(eligiblePts) ==0 %cannot find new eligible spot, just shuffle the rest of points
        p = datasample(availGrid,numel(availGrid),'Replace',false); %Select grids for darts
        tempPts = sGrid(p,:);
        pts = [pts;tempPts];
        break;
    end
        
    %Update empty Grid
    emptyIdx = p(eligiblePts);
    emptyGrid(emptyIdx) = 0;
    
    %Update score pts
    scoreIdx = p(~eligiblePts);
    scoreGrid(scoreIdx) = scoreGrid(scoreIdx) + 1;
    
    %Update emptyGrid if scoreGrid has exceedd k dart throws
    emptyGrid = emptyGrid & (scoreGrid<k);
    
    %Update quantities for next iterations
    nEmptyGrid = sum(emptyGrid);
    pts = [pts;tempPts];
    ptsCreated = size(pts,1);
    iter = iter+1;
end
[~,~,permidx] = intersect(sGrid,pts,'stable','rows'); 
end