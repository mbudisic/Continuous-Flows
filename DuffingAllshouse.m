classdef DuffingAllshouse < ContinuousFlows.ODEFlow2D
%DUFFINGALLSHOUSE Duffing oscillator used in Allshouse, Thiffeault, Chaos
%(2012)

  properties
    alpha % forcing strength in x coordinate
    omega % forcing angular frequency
    gamma % forcing strength in y coordinate
    delta % damping
  end


  methods

    function obj = DuffingAllshouse( dt, flowp )
    %DUFFINGALLSHOUSE Forced, damped Duffing oscillator, as used in
    %Allshouse, Thiffeault (2012), Chaos
    % DuffingAllshouse( dt, params )
    %
    % dt    time discretization step
    % flowp
    %     -- 1 x 4 vector of coefficients [alpha, omega, gamma, delta]
    %              alpha % forcing strength in x coordinate
    %              omega % forcing angular frequency
    %              gamma % forcing strength in y coordinate
    %              delta % damping
    %     -- 'Chaos2012'   - parameter set [0.1, 1, 0.14, 0.08]

      if nargin < 2
        help ContinuousFlows.DuffingAllshouse.DuffingAllshouse
      end

      obj.Domain = [-3,3; -3,3];
      obj.dt = dt;

      if ischar( flowp )
        switch flowp
          case 'Chaos2012'
            flowp = [0.1, 1, 0.14, 0.08];
          otherwise
            error('Unknown parameter set');
        end
      end

      flowp = num2cell(flowp);
      [obj.alpha, obj.omega, obj.gamma, obj.delta] = deal(flowp{:});

      %% Set up integration parameters
      obj.integrator = @ode45;
      obj.intprops = odeset;
      obj.intprops = odeset(obj.intprops, 'Vectorized','on');
    end


    function [ f ] = vf( obj, t, x )
    % VF Compute vector field along a trajectory
    % a single trajectory given by (t, x)
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

      X = x(1,:);
      Y = x(2,:);

      f(1,:) = Y + obj.alpha*cos(obj.omega*t);
      f(2,:) = X-X.^3+obj.gamma*cos(obj.omega*t) - obj.delta*Y;

    end % vf


    function [ J ] = jacobian( obj, t, x )
    % JACOBIAN Compute Jacobian of the vector field along
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
    %     - of the vector field at [ t(i), x(i,:) ] point

      X = x(1,:);
      Y = x(2,:);
      L = numel(t);

      J = nan(2,2,L);

      J(1,1,:) = 0;
      J(1,2,:) = 1;
      J(2,1,:) = 1-3*X.^2;
      J(2,2,:) = -obj.delta;

    end % jacobian

  end % methods

end % classdef
