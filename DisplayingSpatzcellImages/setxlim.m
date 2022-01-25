<<<<<<< HEAD
function maxx = setxlim(varargin) 
% find the maxx so that fraction p of ya is within the range [0 maxx]
% usually for ease of showing a histogram where most the data is 
% xa is usually a bin
% MW, Dec 16, 2016

if numel(varargin) == 3 
    xa = varargin{1} ; 
    ya = varargin{2} ; 
    p = varargin{3} ; 
    i1 = 'prob' ; % a distribution/probability plot
elseif numel(varargin) == 2
    xa = varargin{1} ; 
    p = varargin{2} ; 
    i1 = 'raw' ; % just a raw/scatter plot
end

switch i1
    
    case 'prob'

        p_cov = zeros(size(xa)) ; 
        for i1 = 1:numel(xa)
            p_cov(i1) = sum(ya(xa<=xa(i1))) ./ sum(ya) ; 
        end
        dp = (p_cov - p) ; 
        id = find(dp>0) ; 
        ids = sort(id) ;
        if numel(ids)>1
            maxx = xa(ids(2)) ; 
        elseif numel(ids) == 1
            maxx = xa(ids(1)) ; 
        elseif numel(ids) == 0
            maxx = xa(end);
        end
    
    case 'raw'
        
        xa = sort(xa) ; 
        p_cov = zeros(size(xa)) ; 
        for i1 = 1:numel(xa)
            p_cov(i1) = length(find(xa<=xa(i1))) ./ length(xa) ; 
        end
        dp = (p_cov - p) ; 
        id = find(dp>0) ; 
        ids = sort(id) ;
        if numel(ids)>1
            maxx = xa(ids(2)) ; 
        elseif numel(ids) == 1
            maxx = xa(ids(1)) ; 
        elseif numel(ids) == 0
            maxx = xa(end);
        end
        
=======
function maxx = setxlim(varargin) 
% find the maxx so that fraction p of ya is within the range [0 maxx]
% usually for ease of showing a histogram where most the data is 
% xa is usually a bin
% MW, Dec 16, 2016

if numel(varargin) == 3 
    xa = varargin{1} ; 
    ya = varargin{2} ; 
    p = varargin{3} ; 
    i1 = 'prob' ; % a distribution/probability plot
elseif numel(varargin) == 2
    xa = varargin{1} ; 
    p = varargin{2} ; 
    i1 = 'raw' ; % just a raw/scatter plot
end

switch i1
    
    case 'prob'

        p_cov = zeros(size(xa)) ; 
        for i1 = 1:numel(xa)
            p_cov(i1) = sum(ya(xa<=xa(i1))) ./ sum(ya) ; 
        end
        dp = (p_cov - p) ; 
        id = find(dp>0) ; 
        ids = sort(id) ;
        if numel(ids)>1
            maxx = xa(ids(2)) ; 
        elseif numel(ids) == 1
            maxx = xa(ids(1)) ; 
        elseif numel(ids) == 0
            maxx = xa(end);
        end
    
    case 'raw'
        
        xa = sort(xa) ; 
        p_cov = zeros(size(xa)) ; 
        for i1 = 1:numel(xa)
            p_cov(i1) = length(find(xa<=xa(i1))) ./ length(xa) ; 
        end
        dp = (p_cov - p) ; 
        id = find(dp>0) ; 
        ids = sort(id) ;
        if numel(ids)>1
            maxx = xa(ids(2)) ; 
        elseif numel(ids) == 1
            maxx = xa(ids(1)) ; 
        elseif numel(ids) == 0
            maxx = xa(end);
        end
        
>>>>>>> 88dd4f4ac7eb6f5ae758552b04c47df241263be0
end