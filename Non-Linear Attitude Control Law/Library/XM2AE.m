function AE = XM2AE(XM, e)
% Description:
% Computation of Eccentric Anomaly from Mean Anomaly.
% It is an iterative process -> we must apply Newton's Method

% Note: all angles are considered in [rad]

f = @(AE) AE - e*sin(AE) - XM;
df = @(AE) 1 - e*cos(AE);

i=1;                % starting iteration counter
max_iter = 50;      % set maximum n° of iterations
tolx = 1e-7;        % set tolerance for x

x = zeros(max_iter, 1);     % vector of solutions for each iteration
x(1) = XM;          % appropriate initial guess if e < 0.7

while i < max_iter
   x(i+1) = x(i) - (f(x(i)) / df((x(i))));      % evaluate new value
   eps = abs(x(i+1) - x(i));            % compute abs error

   if eps < tolx            % break cycle if tolerance is satisfied
       AE = x(i+1);
       break
    
   elseif i == max_iter-1   % if we are exceeding max n° iterations
       AE = x(i+1);         

   end
   
   i = i + 1;

end

end