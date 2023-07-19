function [JPdot] = buildJPdot(md, xi, xidot)

% Description: this function builds the derivative in time of JP which is
% equal to the one of JDP since the latter is the only time-varying term in
% the former.

JPdot = md * xidot * [2*xi, 0, 0;...
                      0, 0, 0;...
                      0, 0, 2*xi];

end