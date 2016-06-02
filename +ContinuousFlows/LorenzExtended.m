classdef LorenzExtended < ContinuousFlows.AbstractODEFlow
%LORENZEXTENDED Lorenz system, as extended by Lyubimov, Zaks (Phys D, 1983)
%
% x' = S (y - x) + S D y(z - R)
% y' = Rx - y - xz
% z' = xy - bz + ax
%

  properties
    % parameters for the flow
    S      % Prandtl number (Lorenz)
    R      % Rayleigh  (Lorenz)
    B      % geometric parameter (Lorenz)
    D      % vibrational parameter
    A      % symmetry-breaking parameter
  end

  methods

    function obj = LorenzExtended( dt, params )
    %ABCFLOW Construct the extended Lorenz system object.
    % LorenzExtended( dt, params )
    % Params can be:
    %
    % -- 1 x 5 vector of coefficients [S,R,B,D,A]
    % -- 'lorenz' - butterfly attractor [10, 28, 8/3, 0, 0]
    % -- 'pikovskyA'- singular-continuous spectrum (Pikovsky, 1994)
    % -- 'pikovskyB'- singular-continuous spectrum (Pikovsky, 1994)
    %
    % where:
    % S      % Prandtl number (Lorenz)
    % R      % Rayleigh  (Lorenz)
    % B      % geometric parameter (Lorenz)
    % D      % vibrational parameter
    % A      % symmetry-breaking parameter
    %
    % dt = 0.01 is a good choice.

      obj.dt = dt;
      obj.Domain = [-15,15; -15, 15; -5,30 ];

      if ischar( params )
        switch params
          case 'lorenz'
            params = [10, 28, 8/3, 0, 0];
          case 'pikovskyA'
            params = [10, 15.8237366, 8/3, 0.0526349, 0];
          case 'pikovskyB'
            params = [10, 14.1487968, 8/3, 0.05433476, -0.56112733];

          otherwise
            error('Unknown parameter set');
        end
      end

      params = num2cell(params);
      [obj.S, obj.R, obj.B, obj.D, obj.A] = deal(params{:});

      %% Set up integration parameters
      obj.integrator = @ode45;
      obj.intprops = odeset;
      obj.intprops = odeset(obj.intprops, 'Vectorized', 'on');
      obj.intprops = odeset(obj.intprops, 'MaxStep', dt);
      obj.intprops = odeset(obj.intprops, 'Jacobian', @(t,x)obj.jacobian(t,x) );
      % obj.intprops = odeset(obj.intprops, 'Stats','on' );

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
    %      x = 2*pi*x;


    z = x(3,:);
    y = x(2,:);
    x = x(1,:);

    f(1,:) = obj.S*(y-x) + obj.S*obj.D*y.*(z-obj.R);
    f(2,:) = obj.R*x - y - x.*z;
    f(3,:) = x.*y-obj.B*z+obj.A*x;

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

      J = nan(3,3,L);

      z = x(3,:);
      y = x(2,:);
      x = x(1,:);

      J(1,1,:) = -obj.S;
      J(1,2,:) = obj.S*( 1 + obj.D*(z-obj.R) );
      J(1,3,:) = obj.S*obj.D*y;

      J(2,1,:) = obj.R - z;
      J(2,2,:) = -1;
      J(2,3,:) = -x;

      J(3,1,:) = y + obj.A;
      J(3,2,:) = x;
      J(3,3,:) = -obj.B;

    end

  end

end
