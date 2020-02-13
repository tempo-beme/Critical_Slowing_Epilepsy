function [width, height] = fwhm(x,y, half)
    
if(half)
    half_height = 0.5;
    %intersecting line
    X = [x(1) x(end)];
    Y = [half_height half_height];
    P = InterX([x(:)'; y(:)'], [X(:)'; Y(:)']);
    
    if(isempty(P))
        width = nan;
        height = 0;
        return;
    end
    b = P(1,1);
    width = b;
    height = y(1);
    
     
else
    [mxy idy1] = max(y);
    [mny idy2] = min(y);
    half_height = 0.5*(mxy - mny) + mny;
    
    %intersecting line
    X = [x(1) x(end)];
    Y = [half_height half_height];
    P = InterX([x(:)'; y(:)'], [X(:)'; Y(:)']);
    height = mxy;
    try
        a = find(P(1,:)<x(idy1));a = a(end);
        b = find(P(1,:)>x(idy1));b = b(1);
    catch
        width = nan;
        height = nan;
        return;
    end
    
    
    if(size(P,2) == 1)
        width = nan;
    else
        width = P(1, b) - P(1,a);
    end
     
end

end