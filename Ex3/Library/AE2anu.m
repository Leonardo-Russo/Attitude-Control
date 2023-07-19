function anu = AE2anu(AE, e)
% Description:
% Computation of True Anomaly from knowledge of Eccentric Anomaly and
% Eccentricity.

% Note: all angles are considered in [rad]

laps = 0;
stop = 0;
anu = 0;

while stop == 0

    if AE > 2*pi
        AE = AE - 2*pi;
        laps = laps+1;

    elseif AE < 0
        AE = AE + 2*pi;
        laps = laps-1;

    else 
        stop = 1;
    
    end
end

if AE > 0 && AE < pi
    anu = 2*atan(sqrt((1+e)/(1-e)) * tan(AE/2));

elseif AE == 0 || AE == pi || AE == 2*pi
    anu = AE;

elseif AE > pi && AE < 2*pi
    anu = 2*pi + 2*atan(sqrt((1+e)/(1-e)) * tan(AE/2));

end

end   