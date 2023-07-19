function AE = anu2AE(anu, e)
% Description: 
% Conversion from True Anomaly to Eccentric Anomaly

% We must adjust the values of the anomalies in the case that the value is
% not inside the [0, 2pi] interval.


laps = 0;
stop = 0;
AE = 0;

while stop == 0

    if anu > 2*pi
        anu = anu - 2*pi;
        laps = laps+1;

    elseif anu < 0
        anu = anu + 2*pi;
        laps = laps-1;

    else 
        stop = 1;
    
    end
end

if anu < pi && anu > 0
    AE = 2 * atan(sqrt((1-e)/(1+e)) * tan(anu/2));

elseif anu == pi || anu == 0 || anu == 2*pi 
    AE = anu;

elseif anu > pi && anu < 2*pi
    AE = 2*pi + (2 * atan(sqrt((1-e)/(1+e)) * tan(anu/2)));

end

AE = AE + laps*2*pi;

end

