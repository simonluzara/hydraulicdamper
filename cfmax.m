function Cfmax=cfmax(D_ori,L_ori)

if ~isempty(D_ori)
    Cfmax = zeros(1, length(D_ori));
    ratios = L_ori ./ D_ori;  % Compute ratios once

    % Conditions for flow coefficients using vectorized operations
    Cfmax(ratios > 0 & ratios <= 1) = 0.61 + 0.16 * ratios(ratios > 0 & ratios <= 1).^2;
    Cfmax(ratios > 1 & ratios < 2) = 0.73 + 0.04 * ratios(ratios > 1 & ratios < 2);
    Cfmax(ratios >= 2 & ratios < 10) = 0.827 - 0.0085 * ratios(ratios >= 2 & ratios < 10);
else
    Cfmax = [];
end
end
