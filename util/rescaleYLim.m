function rescaleYLim(h, bounds)
%RESCALEYLIM Rescale axes YLim property for set bounds
%   RESCALEYLIM(h, bound) rescales the axes with handle <h> for the bounds
%   set in <bounds>.
%
%   Arguments:
%   - h : the axes handle
%   - bounds : a 2 element vector containing the lower bound LB and upper
%   bound UB. bounds = [LB,UB]. For no lower bound set LB to -Inf and if
%   there is no upper bound set UB to Inf.
b = get(h, 'YLim');
if(bounds(1) ~= -inf)
    b(1) = bounds(1);
end
if(bounds(2) ~= inf)
    b(2) = bounds(2);
end
set(h, 'YLim', b); 