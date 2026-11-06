function [v_disp, v_vel, v_acc, z, dot_z] = ode45u(dxdt,time,x0,u,params,tolerance,display)
% [time, x, x_dot, y] = ode45u( dxdt, time, x0, u, params, tolerance, display ) 
%
% Solve a system of nonhomogeneous ordinary differential equations
% using the embedded Runge-Kutta formulas of J.R. Cash and A.H. Karp, and
% using linear interpolation of external forcing within time-steps.  
%
%  Input Variable      Description
%  --------------      -----------
%   dxdt          : a function of the form [x_dot,y] = dxdt(t,x,u,params)
%                   which provides the state derivative given the time, t, 
%                   the state, x, external forcing, u, and parameters
%                   a vector system outputs, y, may also be returned.
%   time          : 1-by-p vector of time values, uniformly increasing
%                   (or decreasing), at which the solution is computed.
%   x0            : n-by-1 vector of initial state values, at time(1).
%   u             : m-by-p matrix of system forcing data for the ode, 
%                   sampled at points in the 1-by-p vector time.
%   params        : optional parameters used in the function dxdt
%   tolerance     : the desired tolerance constant (default 0.001)
%                   may be a vector of the length of x0, or a scalar
%   display       : 3= display lots of results; 2= display less; 1= even less 
%
% Output Variable        Description
% ---------------        -----------
%   time          : is returned, un-changed
%   x             : the solution to the differential equation at times
%   x_dot         : the derivative to the differential equation at times
%   y             : the output of the differential equation at times
%

 if nargin < 3
    help ode45u
    disp(' *** Not enough input arguments. ***');
    disp(' *** usage: ode45u( dxdt, time, x0, u, params, tolerance, display )');
    return
 end

 points  = length(time);   % number of points in the time history
 n       = length(x0(:));  % number of elements in the state vector

 if nargin < 4, u = zeros(1,points);   end
 if nargin < 5, params    = 0;         end
 if nargin < 6, tolerance = 1e-3;      end
 if nargin < 7, display   = 0;         end
 if nargin > 5, tolerance = max(tolerance,1e-12);  end
 if nargin > 6, display   = ceil(min(display,3));  end


 [dxdt0,y0]  = feval( dxdt, time(1), x0, u(:,1), params );  % init. output

 m = length(y0(:));

 [ru,cu] = size(u);    

 MaxSteps    = 1000;        % Maximum number of interior steps
 MaxNumSteps = 0;           % Maximum number of sub-steps taken in an interval
 MaxError    = 0;           % Maximum sub-step truncation error 

% Cash-Karp coefficients for embedded Runge-Kutta steps
  a1=    0;
  a2=   1/5;      a3=  3/10;   a4=  3/5;      a5= 1;            a6=  7/8;
 b21=   1/5;
 b31=   3/40;    b32=  9/40;
 b41=   3/10;    b42= -9/10;  b43=  6/5;
 b51= -11/54;    b52=  5/2;   b53=-70/27;    b54=   35/27;
 b61=1631/55296; b62=175/512; b63=575/13824; b64=44275/110592; b65=253/4096;

 x      = zeros(n,points);  % memory allocation for the solution states
 x_dot  = zeros(n,points);  % memory allocation for the solution rates
 y      = zeros(m,points);  % memory allocation for the solution output
 x(:,1) = x0(:);            % the initial conditions
 y(:,1) = y0(:);            % the initial conditions
 XOLD   = x0;               % state vector at the start of each interior step
 fevals = 0;                % running sum of function evaluations
 SingleStep = 0;            % 1: run as a single-step solver; 0: don't

 if cu < points
    u(:,cu+1:points) = 0;   % pad u with zeros if it is not long enough
 end

 for p = 1:points-1   % ====== loop over all points in the time series

     t0 = time(p);                 % time at start of the full step
     dt = time(p+1) - t0;          % size of the full step
     T0 = t0;                      % time at start of an interior sub-step
     DT = dt;                      % size of the interior sub-steps
     
     NumSteps = 1;                 % number of interior steps 
     step     = 0;                 % the current interior step number

     u0   = u(:,p);                % forcing at start of full time step
     dudt = (u(:,p+1)-u0)/dt;      % change in forcing over time step

     while step < NumSteps   % ------ loop to adjust the integration time step 

         [k1,y1] = feval( dxdt, T0 + DT*a1, XOLD, u0+dudt*(T0+DT*a1-t0), params );
         [k2,y2] = feval( dxdt, T0 + DT*a2, XOLD + k1*DT*b21, u0+dudt*(T0+DT*a2-t0), params );
         [k3,y3] = feval( dxdt, T0 + DT*a3, XOLD + k1*DT*b31 + k2*DT*b32, u0+dudt*(T0+DT*a3-t0), params );
         [k4,y4] = feval( dxdt, T0 + DT*a4, XOLD + k1*DT*b41 + k2*DT*b42 + k3*DT*b43, u0+dudt*(T0+DT*a4-t0), params );
         [k5,y5] = feval( dxdt, T0 + DT*a5, XOLD + k1*DT*b51 + k2*DT*b52 + k3*DT*b53 + k4*DT*b54, u0+dudt*(T0+DT*a5-t0), params  );
         [k6,y6] = feval( dxdt, T0 + DT*a6, XOLD + k1*DT*b61 + k2*DT*b62 + k3*DT*b63 + k4*DT*b64 + k5*DT*b65, u0+dudt*(T0+DT*a6-t0), params );

         % 5th order predictor
         x5 = XOLD + ( k1*37/378 + k3*250/621 + k4*125/594 + k6*512/1771 )*DT;
         % 4th order predictor
         x4 = XOLD + ( k1*2825/27648 + k3*18575/48384 + k4*13525/55296 + k5*277/14336 + k6*1/4 )*DT;

         fevals = fevals + 6;

         % evaulate the truncation error at the start of the full step
         if step == 0
            for i=1:n
               if x5(i) ~= 0
                  TruncationError(i) = abs ( x5(i) - x4(i) )/abs(x5(i)); 
               else
                  TruncationError(i) = 0.0; % avoid Inf values
               end
            end
            Converged = all ( TruncationError <= tolerance );
         end
         if ~all(isfinite(TruncationError)), break; end  % floating point error

         if Converged && NumSteps == 1       % a single full dt step is ok
            XOLD = x5;
            break

         elseif ~Converged && step == 0      % look for a good sub-step size
            % Increase NumSteps using eq'n 16.2.10 of Numerical Recipes in C ,
            NS1 = NumSteps / min( (tolerance./TruncationError).^0.25 );
            % ... or increase NumSteps by 10 percent, ...
            NS2 = 1.1 * NumSteps; 
            % ... whichever is greater ...
            NumSteps = ceil( max( NS1 , NS2 ) );  
            DT = dt / NumSteps;              % smaller interior sub-step size

         elseif Converged && step < NumSteps  % the step size is ok, continue
            step = step + 1;                 % increment the sub-step number
            T0 = t0 + step*DT;           % starting time for the next sub-step
            XOLD = x5;                   % state at start of next sub-step

         else
            disp(' uh-oh ... this should never happen ... check ode45u.m ...'); 
            return
         end

         if display == 3  &&  step == 1       % display intermediate information
            fprintf('point   NS   s   error      time     fevals\n');
         end
         if display == 3  &&  step > 0        % display intermediate information
            fprintf('%5d  %3d %3d %9.2e %10.6f   %4d \n', ...
                p, NumSteps, step, max(TruncationError), T0, fevals);
         end

         if NumSteps > MaxSteps              % too many interior steps
            fprintf('ode45u: The required interior step size is too small\n');
            fprintf('time= %10.4f NumSteps= %4d  DT= %9.2e  Error= %9.2e\n',...
                     p*dt, NumSteps, DT, max(TruncationError) );
	    if ~isreal(max(TruncationError))
		return;
	    end
            if tolerance < 1.00
               tolerance = tolerance*1.5;
               fprintf('Increasing the desired tolerance and trying again');
               fprintf(' ... tolerance = %10.2e\n', tolerance);
            elseif MaxSteps < 5e3
               MaxSteps = MaxSteps*1.1;
               fprintf('Increasing maximum allowable steps and trying again');
               fprintf(' ... MaxSteps = %9d\n', MaxSteps);
            else
               fprintf('continuing  tolerance = %10.2e   MaxSteps = %9d\n',...
                                                      tolerance, MaxSteps );
            end
            NumSteps = MaxSteps;
         end

     end % ------ end of loop to adjust the integration time step 
     
     MaxError    = max(max(TruncationError),MaxError);
     if NumSteps > MaxNumSteps
        MaxNumSteps = NumSteps;
        pMax = p;
     end

     if display == 2 && NumSteps > 1          % display intermediate information
        fprintf('p = %6d  NumSteps = %3d  Error = %9.2e\n', ...
                               p, NumSteps, max(TruncationError) );
     end

     x(:,p+1) = XOLD;                          % save the solution state
     x_dot(:,p) = k1(:);                       % save the solution rate
     y(:,p) = y1(:);                           % save the solution output
     if ( norm(k1)*dt > 1e6*norm(XOLD) ) break; end  % break if unstable

 end % ====== end of loop over all points in the time series

 x_dot(:,p+1) = k1(:);                         % save the final solution rate
 y(:,p+1) = y1(:);                             % save the final solution output

 if display > 0         % display overall solution statistics
    fprintf('p = %4d   MaxNumSteps = %8d   MaxError = %10.2e  fevals = %d\n',...
                                    pMax, MaxNumSteps, MaxError, fevals );
                                
 end
 v_disp = x(1,:); v_vel = x_dot(1,:); v_acc = x_dot(2,:);
 z = y; dot_z = x_dot(3,:);
 end
%  References
%
%  Cash, J.R., and Karp, A.H. 1990, 
%  A Variable Order Runge-Kutta Method for Initial-Value Problems with 
%  Rapidly Varying Right-Hand Sides,
%  ACM Transactions on Mathematical Software, vol. 16, pp. 201-222.
%  http://www.duke.edu/~hpgavin/ce283/Cash-90.pdf
%
%  Press, W.H., Teukolsky, S.A., Vetterling, W.T., and Flannery, B.P.,
%  Numerical Recipes in C, Cambridge Univ Press, 1992, (ISBN 0-521-43108-5)
%  Section 16.1-16.2
%  http://www.nrbook.com/

% ------------------------------------------------------------------- ODE45U
