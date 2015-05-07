%%HACKBORNROTOSC
% Hackborn Rotor-Oscillator flow -- Hackborn et al., JFM, (1997)
%
% Flow domain is [-1, 1] x [-ymax, ymax]. In theory, ymax -> infinity. In
% practice, ymax of 4 and more is sufficient.
%
% Part of the velocity field is evaluated at a grid and then interpolated.


classdef HackbornRotOsc < ContinuousFlows.ODEFlow

  properties
    %% flow properties
    epsilon % strength of wall oscillation
    lambda  % frequency of wall oscillation
    c       % rotor location (between -1 and 1)

    %% discretization properties
    ymax    % length of the channel (larger than 1)

    %% quadrature
    quadk   % coordinate points
    quadw   % weights
  end


  methods

    function obj = HackbornRotOsc( dt, flowp, ymax )
    %HACKBORNROTOSC Construct a Hackborn Rotor-Oscillator flow
    % HackbornRotOsc( dt, params )
    %
    % dt    time discretization step
    % flowp is a 1 x 3 vector of flow parameters [epsilon, lambda, c]
    % ymax  is a scalar determining the height of the domain

      obj.dt = dt;

      flowp = num2cell(flowp);
      [obj.epsilon, obj.lambda, obj.c] = deal(flowp{:});
      obj.ymax = ymax;

      %% Set up integration parameters
      obj.integrator = @ode45;
      obj.intprops = odeset;
      obj.intprops = odeset(obj.intprops, 'Vectorized','off');

      %% Generate Gauss-Legendre points on 0-1 interval
      N = 100;
      [Quadk,Quadw] = lgwt(N, 0, 50);

      obj.quadk = Quadk(:).';
      obj.quadw = Quadw(:).';

    end

    function [ f ] = vf( obj, t, x )
    % VF Compute the vector field along a single
    % trajectory given by (t, x)
    % [ f ] = vf( obj, t, x )
    %
    % t   - row-vector of times
    % x   - trajectory
    %     - columns correspond to time steps
    %     - rows correspond to states
    %
    % Returns:
    % f   - evaluation of the vector field
    %     - each f(:,i) is a dim x 1 vector field evaluation
    %     - of the vector field at [ t(i), x(i,:) ] point

    % Non-autonomous term
      nonaut = obj.epsilon*(1 + x(1,:)).*cos(obj.lambda*t);

      % Autonomous-term
      dF = obj.DlogF( x );
      dG = obj.DG(x);
      f = [0 -1; 1 0] * (dF+dG);

      % Full vector field (time dep. comes only in y coordinate)
      f(2,:) = f(2,:) + nonaut;

    end

    function [dG] = DG(obj, x)
    %DG Indefinite-integral term in the autonomous stream function

      [X,K] = meshgrid( x(1,:), obj.quadk );
      [Y,~] = meshgrid( x(2,:), obj.quadk );

      C = obj.c;

      g = (2*cosh(K*C)./(sinh(2*K) + 2*K)) .* ...
          ( tanh(K) .* cosh(K.*X) - X.* sinh(K.*X) ) + ...
          (2*sinh(K*C)./(sinh(2*K) - 2*K)) .* ...
          ( coth(K) .* sinh(K.*X) - X.* cosh(K.*X) ) ;

      gx = (2*cosh(K*C)./(sinh(2*K) + 2*K)) .* ...
           ( K.*tanh(K) .* sinh(K.*X) - K.*X.* cosh(K.*X) - sinh(K.*X) ) + ...
           (2*sinh(K*C)./(sinh(2*K) - 2*K)) .* ...
           ( K.*coth(K) .* cosh(K.*X) - K.*X.* sinh(K.*X) - cosh(K.*X) );

      Gx =  obj.quadw * (gx.*    cos(K.*Y)); % weighted sum as an inner product
                                           % with weight row-vector
      Gy = -obj.quadw * (g .* K.*sin(K.*Y));

      dG = [Gx; Gy];
      assert( all(size(dG) == size(x)) );

    end

    function [dF] = DlogF( obj, x )
    %DLOGF Closed-form term in the autonomous stream function

    % The term in question is 0.5 * log(F)
    % where F can be written as F = (A - B) / (A+B)
    % which in turn means that
    %
    % d log(F)/2 = (A * dA - B * dB)/(A^2 - B^2)

    C = obj.c;

    y = x(2,:);
    x = x(1,:);

    Den = ( cos(pi.*(C-x)/2)-cosh(pi.*y/2) ) .* ( cos(pi.*(C+x)/2)+cosh(pi.*y/2) );

    dFx =  (pi/4).*( sin(C.*pi)-2.*cos(C.*pi/2).*cosh(pi.*y/2).*sin(pi.*x/2) )./Den;
    dFy = -(pi/2).*( cos(C.*pi/2).*cos(pi.*x/2).*sinh(pi.*y/2) )./Den;

    dF = [dFx; dFy];

    end

    function [varargout] = quiver( obj, t, N )
    %QUIVER Vector field of the flow.
    %
    % Produce vector field of the flow in 2*N x N points
    % at time t.
    %
    % QUIVER(obj, t, N)
    %   Plots the vector field at time t on grid of N x N points.
    % h = QUIVER(obj, t, N)
    %   As above, and returns graphics handle of the quiver object.
    % [X,Y,U,V] = QUIVER(obj, t, N)
    %   Returns spatial points and components of the vector field.

      [X,Y] = meshgrid(linspace(-1,1,N), ...
                       linspace(-obj.ymax,obj.ymax,N) );
      f = obj.vf(t, [X(:),Y(:)].');
      U = reshape(f(1,:), size(X));
      V = reshape(f(2,:), size(Y));

      if nargout > 1
        varargout = {X,Y,U,V};
      else
        h = quiver(X,Y,U,V);
        if nargout > 0
          varargout = h;
        end
      end

    end


    function jacobian( obj, t, x )
      pass
    end

  end

end
