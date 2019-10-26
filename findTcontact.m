function [tTakeoff, idxTakeoff] = findTcontact(t,Y,T,th,varargin)
    if nargin == 4
        edge = 'up';
    else
        edge = varargin{1};
    end
    Yshift = [Y(2:end) Y(end)];
    if strcmpi(edge,'up')
        idxTakeoff = find(t > T & Y < th & Yshift >= th, 1,'first');
    else
        idxTakeoff = find(t > T & Y > th & Yshift <= th, 1,'first');
    end
    tTakeoff = t(idxTakeoff);
    
end