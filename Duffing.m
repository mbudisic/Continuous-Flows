classdef Duffing < ContinuousFlows.Hamiltonian2DFlow
%DUFFING Unforced, undamped Duffing oscillator.
%
% Energy function is alpha*x^4/4 + beta*x^2/2 + y^2/2

  properties
    %% flow properties
    alpha % nonlinear spring stiffness
    beta  % linear spring stiffness
    epsilon % strength of force
    omega % frequency of force

  end


  methods

    function obj = Duffing( dt, flowp )
    %DUFFING Construct a 2d unforced undamped Duffing oscillator
    % Duffing( dt, alpha, beta )
    %
    % dt    time discretization step
    % flowp
    %     -- 1 x 3 vector of coefficients [alpha, beta]
    %                  alpha   -- nonlinear spring constant
    %                  beta    -- linear spring constant
    %     -- 'doublewell'   - parameter set [1, -1, 0, 0]

      if nargin < 3
        help ContinuousFlows.Duffing.Duffing
      end

      obj.Domain = [-2,2; -2,2];
      obj.dt = dt;

      if ischar( flowp )
        switch flowp
          case 'doublewell'
            flowp = [1, -1, 0, 0];
          otherwise
            error('Unknown parameter set');
        end
      end

      flowp = num2cell(flowp);
      [obj.alpha, obj.beta, obj.epsilon, obj.omega] = deal(flowp{:});

      %% Set up integration parameters
      obj.integrator = @ode23t;
      obj.intprops = odeset;
      obj.intprops = odeset(obj.intprops, 'Vectorized','on');
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

      Nx = size(x,2);

      X = x(1,:);
      Y = x(2,:);

      if o == 0
        out = obj.beta*0.5*X.^2 + obj.alpha*0.25*X.^4+ 0.5*Y.^2 + X.*obj.epsilon*cos(obj.omega*t);
      elseif o == 1
        out = [obj.beta*X + obj.alpha*X.^3; Y];
      elseif o == 2
        out = [obj.beta + 3*obj.alpha*X.^2; ... % dxxPsi
               zeros(1,Nx); ... % dxyPsi
               ones(1,Nx)]; % dyyPsi
      else
        error('Higher orders not implemented');
      end

    end

  end

end
