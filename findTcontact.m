function [tTakeoff, idxTakeoff] = findTcontact(t,Y,T,th)
    Yshift = [Y(2:end) Y(end)];
    idxTakeoff = find(t > T & Y < th & Yshift >= th, 1,'first');
    tTakeoff = t(idxTakeoff);
end