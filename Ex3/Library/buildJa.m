function Ja = buildJa(Is, a)
% Description: this function creates the A matrix used in the evaluation of
% omegas_dot from the known value of commanded Torque.

Ja = [];

Nw = size(a, 2);

for i = 1 : Nw
    Ja = [Ja, Is*a(:, i)];
end

end