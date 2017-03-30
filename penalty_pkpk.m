function penalty=penalty_pkpk(pos,saturation)
% penalty=penalty_pkpk(pos)
%
% Methods:
% Penalty is a piecewise defined function: low oscs incur no penalty while high oscillations are
% heavily penalised by a power law
%
% TODO: See test_penalty_pkpk.m for behaviour
%
% Modified from AntiCheat.m to penalise position oscillations on detector extending beyond
% detection limits
%
amp_pkpk = max(pos) - min(pos);     % peak-peak amplitude of the captured oscillation

val=amp_pkpk/saturation;        % normalised amplitude

threshold=2/3;
if val<threshold
    penalty=0;
else
    penalty=10*(val-threshold)^4;       % multiplier such that when val==saturation, penalty=0.12
end