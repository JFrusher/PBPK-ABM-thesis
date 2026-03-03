function bsa = calculateBSA(weight_kg, height_cm)
% calculateBSA  Mosteller formula for body surface area (m^2)
%   bsa = sqrt((height_cm * weight_kg) / 3600)
    validateattributes(weight_kg,{'numeric'},{'scalar','positive'});
    validateattributes(height_cm,{'numeric'},{'scalar','positive'});
    bsa = sqrt((height_cm .* weight_kg) / 3600.0);
end