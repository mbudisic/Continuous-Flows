classdef BickleyRossby < ContinuousFlows.Hamiltonian2DFlow
%BICKLEYROSSBY Polar jet with Rossby traveling wave perturbation.
%
% Typical domain size: [0, 20] x [-4,4]
%
% The flow first appears in
% del‐Castillo‐Negrete, Diego, and P. J. Morrison. “Chaotic Transport by Rossby Waves in Shear Flow.” Physics of Fluids A: Fluid Dynamics (1989-1993) 5, no. 4 (April 1, 1993): 948–65. doi:10.1063/1.858639.
%
% and the parameterization used is from
%
% I. Rypina, M.G. Brown, F.J. Beron-Vera, H. Koçak, M.J. Olascoaga, and I.A. Udovydchenkov, J Atmos Sci 64, 3595 (2007).
% "On the Lagrangian Dynamics of Atmospheric Zonal Jets and the
% Permeability of the Stratospheric Polar Vortex"
% doi: 10.1175/JAS4036.1
%
% Stream function:
% $$
% \Psi = c_3 y-u_0 \ell \left(\tanh \left(\frac{y}{\ell
% }\right)-\text{sech}^2\left(\frac{y}{\ell }\right) \sum_{i=1}^3a_i \cos \left(k_i
%   \left(x-c_i t \tau \right)\right)\right)
% $$
%
%
% The stream function is given by
% \begin{align}
%   \Psi(x,y,t) &= c_{3}y - L U_{0} \left\{ \tanh(y/L) -\sech^{2}(y/L) \sum_{i=1}^{3} A_{i} \cos[ k_{i}(x-c_{i} \tau t) ] \right\}  \\
% \shortintertext{with partial derivatives}
% \frac{d\Psi}{dx} &=
% - L U_{0} \sech^{2}(y/L) \sum_{i=1}^{3}
% A_{i} k_{i}\sin[ k_{i}(x-c_{i} \tau t) ] \\
% \frac{d\Psi}{dy} &=
% c_{3} - U_{0}\sech^{2}(y/L)\left\{ 1 + 2 \tanh(y/L)
% \sum_{i=1}^{3} A_{i} \cos[ k_{i}(x-c_{i} \tau t) ]
% \right\}.
% \end{align}
% Vector field is given by
% \begin{equation}
%   f(x,y,t) =
%   \begin{bmatrix*}[r]
%     -\frac{d\Psi}{dy}(x,y,t) \\
%     \frac{d\Psi}{dx}(x,y,t)
%   \end{bmatrix*}.
% \end{equation}
% Constants used in simulation are
% \begin{align}
%   \tau &= 24 \cdot 3600 \cdot 10^{-6} \\
%   U_{0} &= 62.66, \quad
%   r_{e} = 6.371, \quad
%   L = 1.770  \\
%   k_{n} &= \frac{2n}{r_{e}}\\
%   c_{1} &= 0, \quad
%   c_{2} = 0.461 U_{0}, \quad
%   c_{3} = 0.205 U_{0} \\
%   A_{1} &= 0, \quad
%   A_{2} = 0.1, \quad
%   A _{3} = 0.3.
% \end{align}
% Values were taken from Rypina et al., (2007) J Atmos Sci \url{http://dx.doi.org/10.1175/JAS4036.1}, with rescaling to units of mega-meters [Mm] in distance and days in time.


  properties (SetAccess = private)

    tscale = 24*3600/1e6;
    U0 = 62.66;

    Re = 6.371; % unit - 10^6 m
    L = 1.770; % unit - 10^6 m

    c1,c2,c3 %
    k1,k2,k3 %
    A1,A2,A3 %


  end

  methods

    function obj = BickleyRossby( dt, params )
    %RYPINAJET Construct the flow with given time step.
    %
    % flow = BickleyRossby;
    % Recommended resolution is dt = 1e-3;

      obj.Domain = [0,20; -4,4];


      obj.k1 = 2/obj.Re;
      obj.k2 = 4/obj.Re;
      obj.k3 = 6/obj.Re;

      obj.A1 = 0.0;
      obj.A2 = 0.1;
      obj.A3 = 0.3;


      if nargin < 2
        params = 'standard';
      end

      validateattributes(dt, {'numeric'},{'positive'});


      if ischar( params )
        switch params
          case 'bickley'
            obj.c1 = 0.0;
            obj.c2 = 0.0;
            obj.c3 = 0.0;
          case 'standard'
            obj.c1 = 0;
            obj.c2 = 0.205*obj.U0;
            obj.c3 = 0.461*obj.U0;

          otherwise
            error('Unknown parameter set');
        end
      end

      if nargin < 1
        obj.dt = 1e-3;
      else
        obj.dt = dt;
      end

      %% Set up integration parameters
      obj.integrator = @ode23t;
      obj.intprops = odeset;
      obj.intprops = odeset(obj.intprops, 'Vectorized', 'on');
      obj.intprops = odeset(obj.intprops, 'Jacobian', @obj.jacobian);
      obj.intprops = odeset(obj.intprops, 'MaxStep', 1e-1);

    end

    function [out] = Psi( obj, t, x, o )
    % PSI Compute the stream function or its derivatives along a
    % trajectory given by (t, x)
    % [ out ] = Psi( obj, t, x, order )
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

      y = x(2,:);
      x = x(1,:);

      A1COSK1 = obj.A1*cos(obj.k1*(x - obj.c1*t*obj.tscale));
      A2COSK2 = obj.A2*cos(obj.k2*(x - obj.c2*t*obj.tscale));
      A3COSK3 = obj.A3*cos(obj.k3*(x - obj.c3*t*obj.tscale));

      A1SINK1 = obj.A1*sin(obj.k1*(x - obj.c1*t*obj.tscale));
      A2SINK2 = obj.A2*sin(obj.k2*(x - obj.c2*t*obj.tscale));
      A3SINK3 = obj.A3*sin(obj.k3*(x - obj.c3*t*obj.tscale));

      TANHY = tanh(y/obj.L);
      SECHY2 = sech(y/obj.L).^2;

      if o == 0
        out = obj.c3*y - obj.L*obj.U0 * ...
              ( TANHY - ...
                SECHY2 .* ...
                ( A1COSK1 + A2COSK2 + A3COSK3 ) ...
                );
      elseif o == 1
        dFx = -obj.L * obj.U0 * SECHY2 .* ...
              ( obj.k1*A1SINK1 + ...
                obj.k2*A2SINK2 + ...
                obj.k3*A3SINK3  );

        dFy = obj.c3 - obj.U0 * SECHY2 .* ...
              (1 + 2*TANHY.* (A1COSK1 + A2COSK2 + A3COSK3) );

        out = [dFx; dFy];
      elseif o == 2
        dFxx = -obj.U0*obj.L*SECHY2.* ...
               ( obj.k1^2*A1COSK1 + ...
                 obj.k2^2*A2COSK2 + ...
                 obj.k3^2*A3COSK3  );
        dFxy = 2 * obj.U0 * TANHY .* SECHY2 .* ...
               ( obj.k1*A1SINK1 + ...
                 obj.k2*A2SINK2 + ...
                 obj.k3*A3SINK3  );
        dFyy = 2*(obj.U0/obj.L) .* ( SECHY2.*TANHY + ...
                                     SECHY2.^2 .* (cosh(2*y/obj.L)-2) .* ...
                                     (A1COSK1 + A2COSK2 + A3COSK3) );

        out = [dFxx; dFxy; dFyy];
      else
        error('Orders > 2 not implemented');
      end

    end
  end
end
