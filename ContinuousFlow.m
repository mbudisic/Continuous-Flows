%CONTINUOUSFLOW
%  Abstract class specifying interface of a "continuous-time dynamical
%  system".
%

classdef (Abstract) ContinuousFlow

  properties
    % trajectory discretization step
    dt

  end

  methods (Abstract)

    [ t, x ] = traj(obj, x0, T, t0)
    % Compute trajectory from t0 -> t0 + T
    % x0  - initial conditions, each row is an i.c.
    % T  - duration of time
    % t0 - initial time
    %
    % Returns:
    % t  - column-vector of time instances
    % x  - set of trajectories
    %      1st ind - time index
    %      2st ind - dimension of state
    %      3rd ind - trajectory

    [ x ] = flow(obj, x0, T, t0)
    % Evaluate flow map from t0 -> t0 + T
    % x0  - initial conditions, each row is an i.c.
    % T   - duration of time
    % t0  - initial time
    %
    % Returns:
    % t   - column-vector of time instances
    % x   - set of points, of the same shape as x0

    [ J ] = jacobian( obj, t, x )
    % Compute Jacobian of the vector field along
    % a single trajectory given by (t, x)
    %
    % t   - column vector of times
    % x   - trajectory
    %     - rows correspond to time steps
    %     - columns correspond to states
    %
    % Returns:
    % J   - Jacobians
    %     - each J(:,:,i) is a dim x dim Jacobian matrix
    %     - of the vector field at [ t(i), x(i,:) ] point

  end


end
