classdef CherryHamiltonianFlow < ContinuousFlows.AbstractHamiltonian2DFlow
%CHERRYHAMILTONIANFLOW Cherry Flow as described in Zaks (2001) for
%Hamiltonian set of parameters (K=1, b free)
%
% Standard domain: [0,2*pi] x [0, 2*pi] (as a 2-torus)
%
% Stream function is
%
% Psi(x,y,t) = -bx + y - cos(y) -sin(x) + sin(x)cos(y)
%
% The flow is steady, but because of the topology of the domain, it can be
% perturbed to interesting behavior (see CherryFlow)
%
% Parameter b can be interpreted as the rotation number of points outside
% the vortex.


  properties
    b        % rotation number of points outside the vortex
  end

  methods

    function obj = CherryHamiltonianFlow( dt, b )
    %CHERRYHAMILTONIANFLOW Construct a CherryHamiltonianFlow object.
    % CherryHamiltonianFlow( dt, b )
    %
    % b - rotation number of trajectories outside the vortex
    %     range abs(b) < 1.56644
    %
    %
    % Stream function is
    % Psi(x,y,t) = A sin( pi * f(x,t) ) sin( pi * y )
    %
    % where f(x,t) = a(t) x^2 + b(t) x
    % a(t) = epsilon * sin( omega t )
    % b(t) = 1 - 2 * epsilon * sin( omega t )
    %


      if nargin < 2
        help ContinuousFlows.CherryHamiltonianFlow.CherryHamiltonianFlow
      end

      obj.dt = dt;
      obj.Domain = [0, 2*pi; 0,2*pi];

      validateattributes(b, {'numeric'}, {'real','scalar'});

      obj.b = b;

      %% Set up integration parameters
      obj.integrator = @ode23t;
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

      if o == 0
        out = -obj.b*X + Y - cos(Y) - sin(X) + sin(X).*cos(Y);
      elseif o == 1
        dXPsi = -obj.b - cos(X) + cos(X).*cos(Y);
        dYPsi = 1 + sin(Y) - sin(X).*sin(Y);
        out = [dXPsi; dYPsi];
      elseif o == 2
        dXXPsi = sin(X) - sin(X).*cos(Y);
        dYYPsi = cos(Y) - sin(X).*cos(Y);
        dXYPsi = -cos(X).*sin(Y);
        out = [dXXPsi; dXYPsi; dYYPsi];
      else
        error('Higher orders not implemented');
      end
    end
  end

end
