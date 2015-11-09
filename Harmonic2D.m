classdef Harmonic2D < ContinuousFlows.Hamiltonian2DFlow
%HARMONIC2D Simple linear harmonic oscillator. Energy function is x^2/2 + y^2/2

  methods

    function obj = Harmonic2D( dt )
    %HARMONIC2D Construct a 2d Harmonic oscillator
    % Harmonic2D( dt )
    %
    % dt    time discretization step

      obj.dt = dt;

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

      if o == 0
        out = 0.5*x(1,:).^2 + 0.5*x(2,:).^2;
      elseif o == 1
        out = x;
      elseif o == 2
        out = [ones(1,Nx); zeros(1,Nx); ones(1,Nx)];
      else
        error('Higher orders not implemented');
      end

    end

  end

end