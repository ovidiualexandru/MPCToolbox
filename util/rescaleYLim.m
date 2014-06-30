function rescaleYLim(h, bounds)
b = get(h, 'YLim');
if(bounds(1) ~= -inf)
    b(1) = bounds(1);
end
if(bounds(2) ~= inf)
    b(2) = bounds(2);
end
set(h, 'YLim', b); 