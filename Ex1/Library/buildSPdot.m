function [SPdot] = buildSPdot(md, xidot)

% Description: this function builds the derivative in time of SP which is
% equal to the one of SDP since the latter is the only time-varying term in
% the former.

SPdot = md * [0, xidot, 0]';

end