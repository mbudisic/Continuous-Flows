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
    %
    % Returns:
    % t  - row-vector of time instances
    % x  - set of trajectories
    %      1st ind - dimension of state
    %      2st ind - time index
    %      3rd ind - trajectory

    % decide if full trajectory is returned or just the last point
      fulltraj = nargout >= 2;

      if nargin < 4
        t0 = 0;
      end

      % initialize output structures
      M = size(x0, 1);
      N = size(x0, 2);
      if fulltraj
        t = ( t0:obj.dt:(t0+T) );
        L = numel(t);
        x = nan( M, L, N );
      else
        M = size(x0, 1);
        N = size(x0, 2);
        x = nan( M, N );
      end

      % timing included
      if nargout >= 3
        timeme = true;
        timings = zeros(1,N);
      else
        timeme = false;
      end

      %%
      % Integrate initial conditions
      for n = 1:N
        fprintf('Computing initial condition %d\n',n)
        disp(x0(:,n));
        if timeme
          tic;
        end
        sol = obj.integrator( @obj.vf, [t0, t0+T], x0(:,n), obj.intprops );
        if timeme
          timings(n) = toc;
        end

        if fulltraj
          % record full trajectory
          x(:,:, n) = deval( sol, t );
        else
          % record just last point
          x(:,n) = sol.y(:,end);
        end
        fprintf(1,'Initial condition %d computed\n', n);
      end

      %%
      % Assign outputs
      varargout{1} = x;
      if fulltraj
        varargout{2} = t;
      end
      if timeme
        varargout{3} = timings;
      end

    end

  end

end
