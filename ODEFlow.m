%CONTINUOUSFLOW
%  Abstract class specifying interface of a "continuous-time dynamical
%  system".
%

classdef (Abstract) ODEFlow < ContinuousFlows.ContinuousFlow

  properties
    integrator
    intprops
  end

  methods

    function [ varargout ] = flow(obj, x0, T, t0)
    % FLOW Compute trajectory from t0 -> t0 + T
    %
    % [ x, t ] = obj.flow(x0, T, t0)
    % x0  - initial conditions, each column is an i.c.
    % T  - duration of time
    % t0 - initial time
    % x = FLOW( x0, T, t0 )
    %     Calculate the values of trajectories at t0+T, for the initial
    %     condition (t0, x0)
    % [ x, t ] = obj.flow(x0, T, t0)
    %     Calculate the full evolution of trajectories the initial
    %     condition (t0, x0) until t0+T, sampled using dt of the object.
    %      1st ind - dimension of state
    %      2st ind - time index
    %      3rd ind - trajectory
    % [ x, t, s ] = obj.flow(x0, T, t0)
    %     As above, but also return the solution (non-interpolated) object.
    %
    % Returns:
    % t  - row-vector of time instances
    % x  - set of trajectories

    % decide if full trajectory is returned or just the last point
      returnFulltraj = nargout >= 2;
      returnSolutions = nargout >= 3;

      if nargin < 4
        t0 = 0;
      end

      % initialize output structures
      M = size(x0, 1);
      N = size(x0, 2);
      t = ( t0:obj.dt:(t0+T) );
      L = numel(t);
      x = nan( M, L, N );
      xf = nan( M, N ); % final point

      % timing included

      %%
      % Integrate initial conditions
      fprintf('Running %d initial conditions\n',N);
      s = cell(N,1);
      parfor n = 1:N
        sol = obj.integrator( @obj.vf, [t0, t0+T], x0(:,n), obj.intprops );

        if returnFulltraj
          % record full trajectory
          x(:,:, n) = deval( sol, t );
        end
        % record just last point
        xf(:,n) = sol.y(:,end);
        fprintf('.');
        if returnSolutions
          s{n} = sol;
        end
      end
      s = [s{:}];
      fprintf('done\n');

      %%
      % Assign outputs
      if returnFulltraj
        varargout{1} = x;
        varargout{2} = t;
      else
        varargout{1} = xf;
      end
      if returnSolutions
        varargout{3} = s;
      end

    end

  end

end
