function [F] = calcPFH(surface1, varargin)
        [tmp Norm2 F] = pclfeature([surface1.X(:),surface1.Y(:),surface1.Z(:)]', 10);
        F = F';
        F  = normalize(F, 'l1', 2);
end
