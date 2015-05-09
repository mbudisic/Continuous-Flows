%CONTINUOUSFLOW
%  Abstract class specifying interface of a "continuous-time dynamical
%  system".
%
%  Leaves open situation where flow map is known explicitly and vector
%  field can be computed from it. For numerical integration of a vector
%  field, see ODEFlows subclass.

classdef (Abstract) ContinuousFlow

  properties
    % trajectory discretization step
    dt

  end

  methods (Abstract)

    [ varargout ] = flow(obj, x0, T, t0)
    %TRAJ Compute trajectory from t0 -> t0 + T
    % [ t, x ] = flow(obj, x0, T, t0)
    %    [ x ] = flow(obj, x0, T, t0)
    %
    % x0  - initial conditions, each column is an i.c.
    % T  - duration of time
    % t0 - initial time
    %
    % If only one output is requested, returns:
    % x   - set of points, of the same shape as x0
    %
    % If two outputs are requested, returns
    % t  - row-vector of time instances
    % x  - set of trajectories
    %      1st ind - dimension of state
    %      2st ind - time index
    %      3rd ind - trajectory

    [ f ] = vf( obj, t, x )
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

    [ J ] = jacobian( obj, t, x )
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

  end


  methods

    function err = testJacobian( obj, t, x, delta )
    %TESTJACOBIAN Compute difference between numeric and analytic
    %             Jacobian matrix.
    %
    % The intended use is to verify correctness of analytic expressions
    % of the Jacobian matrix.
    %
    % err = obj.testJacobian(t,x)
    %
    % Compute the difference between obj.jacobian and second-order central
    % difference of obj.vf at a single space-time point (t,x).
    %
    % err = obj.testJacobian(..., delta)
    %
    % Use spatial step delta in central difference (default: 1e-6)
    %

      if nargin < 4
        delta = 1e-6;
      end

      Nx = size(x,2);
      D  = size(x,1);

      assert( Nx == 1, 'Single point x has to be provided');
      assert( numel(t) == 1, 'Single point t has to be provided');

      %% analytic jacobian
      aJ = obj.jacobian( t, x );

      %% central difference
      stencil = eye(D)*delta;
      xi = [bsxfun(@plus, x, stencil), ...
            bsxfun(@minus, x, stencil) ];
      ti = repmat( t, [1, size(xi,2)] );

      % vector field at stencil points
      v = obj.vf(ti, xi);

      % central difference step
      nJ = ( v(:,1:D) - v(:,(D+1):end) )/(2*delta);

      err = aJ-nJ;

    end

  end


end
