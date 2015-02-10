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


end
