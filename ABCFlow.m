%ABCFLOW
%  Unsteady ABC flow.
%
% dx = ( A + D t sin( pi t ) ) sin(z) + C cos(y)
% dy = B sin(x) + ( A + D t sin( pi t ) ) cos(z)
% dz = C sin(y) + B cos(x)

classdef ABCFlow < dynamics.ContinuousFlow

  properties
    % parameters for the flow
    A
    B
    C
    D
  end

  properties (SetAccess = private )
    integrator
    intprops
  end

  methods

    function obj = ABCFlow( dt, params )
    %ABCFLOW Construct an ABC flow object.
    %
    % Params can be:
    %
    % -- 1 x 4 vector of coefficients [A,B,C,D]
    % -- 'integrable' - integrable parameter set
    % -- 'steady' - chaotic autonomous parameter set
    % -- 'unsteady' chaotic non-autonomous parameter set

      obj.dt = dt;

      if ischar( params )
        switch params
          case 'integrable'
            params = [2, 1, 0, 0];
          case 'steady'
            params = [sqrt(3), sqrt(2), 1, 0];
          case 'unsteady'
            params = [sqrt(3), sqrt(2), 1, 1];
          otherwise
            error('Unknown parameter set');
        end
      end

      params = num2cell(params);
      [obj.A, obj.B, obj.C, obj.D] = deal(params{:});

      %% Set up integration parameters
      obj.integrator = @ode23t;
      obj.intprops = odeset;
      obj.intprops = odeset(obj.intprops, 'Vectorized', 'on');
      obj.intprops = odeset(obj.intprops, 'Jacobian', @(t,x)obj.jacobian(t,x) );
      %      obj.intprops = odeset(obj.intprops, 'Stats','on' );

    end

    function [ f ] = vf( obj, t, x )
    % Compute vector field along
    % a single trajectory given by (t, x)
    %
    % t   - row-vector of times OR a scalar
    % x   - trajectory
    %     - rows correspond to time steps
    %     - columns correspond to states
    %
    % Returns:
    % f   - evaluation of the vector field
    %     - each f(:,i) is a dim x 1 vector field evaluation
    %     - of the vector field at [ t(i), x(i,:) ] point

    % periodicity on the [0,1]^3 cube
      x = 2*pi*x;

      tcoeff = obj.A + obj.D*t.*sin(pi*t);
      f(1,:) = tcoeff .* sin(x(3,:)) + obj.C .* cos(x(2,:));
      f(2,:) = tcoeff .* cos(x(3,:)) + obj.B .* sin(x(1,:));
      f(3,:) = obj.C*sin(x(2,:)) + obj.B*cos(x(1,:));

      % periodicity on the [0,1]^3 cube
      f = f/(2*pi);

    end

    function [ J ] = jacobian( obj, t, x )
    %JACOBIAN Compute Jacobian of the vector field along
    % a single trajectory given by (t, x)
    %
    % [ J ] = obj.jacobian( t, x )
    %
    % t   - row-vector of times OR a scalar
    % x   - trajectory
    %     - columns correspond to time steps
    %     - rows correspond to states
    % Returns:
    % J   - Jacobians
    %     - each J(:,:,i) is a dim x dim Jacobian matrix
    %     - of the vector field at [ t(i), x(i,:) ] point

      assert( numel(size(x)) == 2 );
      L = size(x,2);
      assert( numel(t) == 1 || numel(t) == L, ...
              ['Time is either a scalar or'...
               'has to match number of steps'] );

      tcoeff = obj.A + obj.D*t.*sin(pi*t);

      % scaling the domain does not change the jacobian
      J(1,1,:) =  zeros(1,1,L);
      J(1,2,:) = -obj.C * sin( x(2,:) );
      J(1,3,:) =  tcoeff .* cos( x(3,:) );

      J(2,1,:) = -obj.B * cos( x(1,:) );
      J(2,2,:) =  zeros(1,1,L);
      J(2,3,:) = -tcoeff .* sin( x(3,:) );

      J(3,1,:) = -obj.B * sin( x(1,:) );
      J(3,2,:) =  obj.C * cos( x(2,:) );
      J(3,3,:) =  zeros(1,1,L);

    end

    function [ varargout ] = flow(obj, x0, T, t0)
    % TRAJ Compute trajectory from t0 -> t0 + T
    %
    % [ x, t ] = obj.traj(x0, T, t0)
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
      fulltraj = nargout == 2;

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

      %%
      % Integrate initial conditions
      for n = 1:N
        sol = obj.integrator( @obj.vf, [t0, t0+T], x0(:,n), obj.intprops );
        if fulltraj
          % record full trajectory
          x(:,:, n) = deval( sol, t );
        else
          % record just last point
          x(:,n) = sol.y(:,end);
        end
      end

      %%
      % Assign outputs
      varargout{1} = x;
      if fulltraj
        varargout{2} = t;
      end

    end

  end

end
