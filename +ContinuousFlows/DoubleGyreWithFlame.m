classdef DoubleGyreWithFlame < ContinuousFlows.AbstractHamiltonian2DFlowWithFlame
%DOUBLEGYRE Planar double-gyre system, as used by Shadden (2005)
%
% Standard domain: [0,2] x [0,1]
%
% Stream function is
%
% Psi(x,y,t) = A sin( pi * f(x,t) ) sin( pi * y )
%
% where f(x,t) = a(t) x^2 + b(t) x
% a(t) = epsilon * sin( omega * t )
% b(t) = 1 - 2 epsilon sin( omega t )

  properties
    A        % magnitude of the stream function
    omega    % angular frequency of the perturbation
    epsilon  % magnitude of the perturbation
  end

  methods

    function obj = DoubleGyreWithFlame( dt, varargin )
    %DOUBLEGYREWITHFLAME Construct a Double Gyre object with flame front.
    % DoubleGyre( dt, A, omega, epsilon, v0, delta )
    % Params can be
    %
    % Standard domain: [0,2] x [0,1]
    %
    % Stream function is
    % Psi(x,y,t) = A sin( pi * f(x,t) ) sin( pi * y )
    %
    % where f(x,t) = a(t) x^2 + b(t) x
    % a(t) = epsilon * sin( omega t )
    % b(t) = 1 - 2 * epsilon * sin( omega t )
    %
    %  - v0 is the (unimpeded) speed of the front
    %  - delta is the width of the boundary layer
    %    in which the front speed is tapered to zero


      if nargin < 2
        help ContinuousFlows.DoubleGyreWithFlame.DoubleGyreWithFlame
      end

      obj.dt = dt;
      obj.Domain = [0,2; 0,1; -pi,pi];
      [obj.A, obj.omega, obj.epsilon,obj.v0,obj.delta] = deal(varargin{:});

      %% Set up integration parameters
      obj.integrator = @ode113;
      obj.intprops = odeset;

      obj.intprops = odeset(obj.intprops, 'Vectorized', 'on');
      obj.intprops = odeset(obj.intprops, 'Jacobian', @(t,x)obj.jacobian(t,x) );
      %      obj.intprops = odeset(obj.intprops, 'Stats','on' );

  end



    function [out] = Psi( obj, t, x, o )
    % PSI Compute the stream function or its derivatives along a
    % trajectory given by (t, x)
    % [ f ] = Psi( obj, t, x, order )
    %
    % t   - row-vector of times
    % x   - trajectory
    %     - columns correspond to time steps
    %     - rows correspond to states
    % o   - order of calculation
    %
    % Returns:
    % out   - evaluation of the stream function or its derivatives
    %       - if o == 0, out is 1 x Nx row vector
    %       - if o == 1, out is 2 x Nx row vector; rows are x and y derivatives
    %                    respectively
    %       - if o == 2, out is 3 x Nx row vector; rows are xx, xy, yy
    %                    derivatives respectively

    %      Nx = size(x,2);

      X = x(1,:);
      Y = x(2,:);

      at = obj.epsilon .* sin( obj.omega * t);
      bt = 1 - 2*at;
      F = at .* X .^2 + bt .* X;
      dF = 2*at .*X + bt;
      ddF = 2*at;
      COSF = cos(pi*F);
      SINF = sin(pi*F);
      COSY = cos(pi*Y);
      SINY = sin(pi*Y);

      if o == 0
        out = obj.A.*SINF.*SINY;
      elseif o == 1
        dXPsi = pi*obj.A*SINY.*COSF.*dF;
        dYPsi = pi*obj.A*COSY.*SINF;
        out = [dXPsi; dYPsi];
      elseif o == 2
        dXXPsi = pi*obj.A*SINY.* ...
                 ( ddF .* COSF - pi .* dF .^2 .* SINF );
        dXYPsi = pi^2*obj.A*COSY.*COSF.*dF;
        dYYPsi = -obj.A*pi^2*SINY.*SINF;
        out = [dXXPsi; dXYPsi; dYYPsi];
      else
        error('Higher orders not implemented');
      end
    end
  end

end
