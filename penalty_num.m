function penalty=penalty_num(numcounts,saturation)
% penalty=penalty_num(numcounts,saturation)
%
% Methods:
% penalty_num is a piecewise defined function: high numbers incur no penalty while low counts are
% heavily penalised by a power law
%
% TODO: See test_penalty_num.m for behaviour
%

% amp_pkpk = max(pos) - min(pos);     % peak-peak amplitude of the captured oscillation
num_total=sum(numcounts);

% val=amp_pkpk/saturation;        % normalised amplitude
val=num_total/saturation;        % normalised number

threshold=1;
if val>=threshold
    penalty=0;
else
    penalty=0.5*(threshold-val)^4;       % multiplier such that when val==saturation, penalty=0.12
end