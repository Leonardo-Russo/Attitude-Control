function XM = AE2XM(AE, e)
% Description:
% Conversion from Eccentric Anomaly to Mean Anomaly

XM = AE - e*sin(AE);

end

