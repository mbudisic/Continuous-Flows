classdef PointVortex < ContinuousFlows.AbstractODEFlow
    %POINTVORTEX Point Vortices in Plane moving under Biot-Savart Law
    %
    % dp_i/dt = sum (gamma_j/2/pi) J ( p_i - p_j )/ abs(p_i-p_j)^2
    %
    % where J = [0 -1; 1 0], gamma are signed vorticities
    %
    
    properties
        % vector of signed vorticities
        gamma(1,:) {mustBeNumeric, mustBeNonempty} = nan;
        core(1,1) {mustBePositive, mustBeNonempty} = 1e-6; % radius of core of the vector (no force)
        type (1,1) logical = true; % if true - usual Biot-Savart, if false, replace sum by max.
    end
    
    properties (Dependent)
        N
    end
    
    methods
        
        function obj = PointVortex( dt, gamma )
            obj.dt = dt;
            obj.Domain = [-1,1];
            obj.gamma = gamma;
            
            %% Set up integration parameters
            obj.integrator = @ode113;
            obj.intprops = odeset;
            obj.intprops = odeset(obj.intprops, 'Vectorized', 'off');
            %obj.intprops = odeset(obj.intprops, 'Jacobian', @(t,x)obj.jacobian(t,x) );
        end
        
        function N = get.N(obj)
            N = length(obj.gamma);
        end
        
        function [dr] = vf( obj, t, r )
            dr = obj.biotsavart( r, r );
        end
        
        function [vx, Psix,V] = biotsavart( obj, px, py )
            % BIOTSAVART Use Biot-Savart law to determine velocity at px
            % based on vortices placed at py
            
            % position vectors are now row vectors in a matrix
            rx = reshape( px, [], 2 );
            ry = reshape( py, [], 2 );
            
            % position vectors are now row/column tensors
            % where the coordinate index is the third dimension
            Rx = permute( rx, [1,3,2] );
            Ry = permute( ry, [3,1,2] );
            
            % this allows us to vectorize taking distances between vectors
            D = Rx - Ry; % Nx x Ny x 2 - vectors pointing from y to x
            L = vecnorm(D, 2, 3);
            
            % regularization - velocity is zero inside the vortex (should
            % probably implement something else) - this amounts to
            % "capping" the stream function
            sel = L < obj.core;
            L(sel) = -inf;
            
            % Compute the stream function            
            Psi = -log(L).*obj.gamma/2/pi;
            if ~obj.type
                [Psix, idx] = min(Psi, [], 2);
            else
                Psix = sum(Psi, 2);
                
            end                
            Psix = squeeze( Psix );
                        

            % Now, Biot-Savart produces velocity vectors orthogonal to D
            V = ( cat(3, -D(:,:,2), D(:,:,1) ) )./(L.^2) ...
                .* obj.gamma /2/pi;
            
            if obj.type
                Vx = sum( V, 2 );
            else
                Vx = V(:,idx,:);
            end
            Vx = permute( Vx, [1,3,2] );
            vx = reshape( Vx, [], 1); % make into column
            
            
            
        end
    end
end
