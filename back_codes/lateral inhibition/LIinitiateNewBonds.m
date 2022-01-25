function g = LIinitiateNewBonds(g)
% initiates delta-notch levels for new bonds that may have formed via T1 transitions
nb = length(g.bonds);
nb_prev = length(g.LImodel.bond_notch);
g.LImodel.bond_delta(nb_prev+1:nb) = 0;
g.LImodel.bond_notch(nb_prev+1:nb) = 0;
end