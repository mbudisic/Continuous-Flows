classdef (Abstract) AbstractODEFlow < ContinuousFlows.AbstractContinuousFlow
%ODEFLOW Abstract class specifying interface of a "continuous-time dynamical system".
%

  properties
    integrator
    intprops
  end

  methods

    function [J] = jacobian( obj, t, x, delta)
    % JACOBIAN Compute Jacobian of the velocity field along
    % a single trajectory given by (t, x)
    % [ J ] = jacobian( obj, t, x )
    %
    % t   - row-vector of times
    % x   - trajectory
    %     - columns correspond to time steps
    %     - rows correspond to states
    % Returns:
    % J   - Jacobians
    %     - each J(:,:,i) is a dim x dim Jacobian matrix
    %     - of the velocity field at [ t(i), x(i,:) ] point
    %
    % [ J ] = jacobian( ..., delta )
    %
    % Use delta as the central difference step. If omitted, delta=1e-6.

      if nargin < 4
        delta = 1e-6;
      end

      Nx = size(x,2);
      Dim = size(x,1);
      if isscalar(t)
        t = repmat(t,[1,Nx]);
      end

      J = nan(Dim,Dim,Nx);

      % for each point
      for idx = 1:Nx

        xx = x(:,idx);
        tt = t(idx);

        %% central difference
        stencil = eye(Dim)*delta;
        xi = [bsxfun(@plus, xx, stencil), ...
              bsxfun(@minus, xx, stencil) ];

        ti = repmat( tt, [size(xi,2),1] );
        fi = obj.vf(ti, xi);

        fiD = ( fi( :, 1:Dim) - fi( :, (1:Dim)+Dim) )/(2*delta);

        J(:,:,idx) = fiD.';
      end

      disp('Done')

    end

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
      if ~obj.quiet, fprintf('Running %d initial conditions\n',N); end
      s = cell(N,1);
      parfor n = 1:N
        myic = x0(:,n);
        sol = obj.integrator( @obj.vf, [t0, t0+T], myic, obj.intprops );

        assert(max(sol.x) >= max(t), ...
               sprintf('Simulation for IC %s did not converge.',...
                       mat2str(myic)));

        if returnFulltraj
          % record full trajectory
          x(:,:, n) = deval( sol, t );
        end
        % record just last point
        xf(:,n) = sol.y(:,end);
        if ~obj.quiet
          fprintf('.');
          if mod(n,10) == 0, fprintf('\n'); end
        end
        if returnSolutions
          s{n} = sol;
        end
      end
      s = [s{:}];
      if ~obj.quiet, fprintf('done\n'); end;

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
