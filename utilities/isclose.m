function yesno = isclose(x,y,varargin)
%isclose: when isequal gives 0 but you know the vectors are close within a
%         tolerance
%   x,y       = numbers or vectors
%   tolerance = the difference between numbers should be less than
    if nargin < 3
        tolerance = 10^-10;
    else
        tolerance = varargin{1};
    end
    
    
    yesno = isequal(size(x),size(y)) ...
            && ...
            all(abs(x(:)-y(:)) < tolerance);
end

