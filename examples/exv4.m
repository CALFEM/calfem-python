% CALFEM Vedo Visuaization example (exv4)
% Author: Andreas Åmand

clear
clc

% --- Gauss points & shape functions from soli8e/soli8s ---
% This is used to get stresses at nodes & modal analysis later

g1=0.577350269189626;
gp(:,1)=[-1; 1; 1;-1;-1; 1; 1;-1]*g1;
gp(:,2)=[-1;-1; 1; 1;-1;-1; 1; 1]*g1;
gp(:,3)=[-1;-1;-1;-1; 1; 1; 1; 1]*g1;

xsi=gp(:,1);  eta=gp(:,2); zet=gp(:,3);

%N(:,1)=(1-xsi).*(1-eta).*(1-zet)/8;  N(:,5)=(1-xsi).*(1-eta).*(1+zet)/8;
%N(:,2)=(1+xsi).*(1-eta).*(1-zet)/8;  N(:,6)=(1+xsi).*(1-eta).*(1+zet)/8;
%N(:,3)=(1+xsi).*(1+eta).*(1-zet)/8;  N(:,7)=(1+xsi).*(1+eta).*(1+zet)/8;
%N(:,4)=(1-xsi).*(1+eta).*(1-zet)/8;  N(:,8)=(1-xsi).*(1+eta).*(1+zet)/8;

% ---



d=0.1; % Elements are 100 mm x 100 mm x 100 mm 

% No. elements per direction
nel_x = 50;
nel_y = 4;
nel_z = 3;
nel = nel_x*nel_y*nel_z;

% No. nodes per direction
nnode_x = nel_x+1;
nnode_y = nel_y+1;
nnode_z = nel_z+1;
nnode = nnode_x*nnode_y*nnode_z;

% --- Creates Coord matrix ---

coord = zeros(nnode,3);
row = 1;
for z = 0:nnode_z-1
    for y = 0:nnode_y-1
        for x = 0:nnode_x-1
            coord(row,:) = [x*d,y*d,z*d];
            row = row+1;
        end
    end
end

% --- Creates Dof matrix ---

dof = zeros(nnode,3);
it = 1;
for row = 1:nnode
    for col = 1:3
        dof(row,col) = it;
        it = it + 1;
    end
end
ndof = size(dof,1)*size(dof,2);

% --- Creates Edof and Boundary condition matrices ---
% Boundary conditions: nodes at x = 0m & x = 50m have a displacement of 0

x_step = 1; % Next node in x-direction
y_step = nnode_x; % Next node in y-direction
z_step = (y_step)*nnode_y; % Next node in z-direction

it = 1; % Element number for loops (used as index in edof)
bc_it = 1; % Iteration for bc (used as index in bc)
force_dof_it = 1; % Iteration for point load dofs (used as index in force_dofs)
node = 1; % For keeping track of node

edof = zeros(nel_x*nel_y*nel_z,8*3+1);
bc = zeros(nnode_y*nnode_z*2*3,2);
force_dofs = zeros(25+1,1); % for saving dofs to apply point loads to
for col = 0:nel_z-1 % Loops through z-axis
    node = 1+z_step*col;
    for row = 0:nel_y-1 % Loops through y-axis
        for el = 0:nel_x-1 % Loops through x-axis
            edof(it,1) = it; % Element number, first row in Edof

            % --- First node ---
            edof(it,2:4) = dof(node,:); % Dofs for first element node
            if el == 0  % If element is at x = 0, save bc
                bc(bc_it,1) = dof(node,1);
                bc_it = bc_it + 1;
                bc(bc_it,1) = dof(node,2);
                bc_it = bc_it + 1;
                bc(bc_it,1) = dof(node,3);
                bc_it = bc_it + 1;
            end
            
            % --- Second node ---
            node = node+x_step; % Gets node number
            edof(it,5:7) = dof(node,:); % Gets dofs for node
            if el == nel_x-1 % If element is at x = 5, save bc
                bc(bc_it,1) = dof(node,1);
                bc_it = bc_it + 1;
                bc(bc_it,1) = dof(node,2);
                bc_it = bc_it + 1;
                bc(bc_it,1) = dof(node,3);
                bc_it = bc_it + 1;
            end
            
            % --- Third node ---
            node = node+y_step;
            edof(it,8:10) = dof(node,:);
            if el == nel_x-1
                bc(bc_it,1) = dof(node,1);
                bc_it = bc_it + 1;
                bc(bc_it,1) = dof(node,2);
                bc_it = bc_it + 1;
                bc(bc_it,1) = dof(node,3);
                bc_it = bc_it + 1;
            end

            % If elements at x = 0 and top row, save y-dofs for later
            if (col == 0 && row == 3 && ismember(el,(13:38)) == 1)
                force_dofs(force_dof_it) = dof(node,2);
                force_dof_it = force_dof_it + 1;
            end
            
            % --- Fourth node ---
            node = node-x_step;
            edof(it,11:13) = dof(node,:);
            if el == 0
                bc(bc_it,1) = dof(node,1);
                bc_it = bc_it + 1;
                bc(bc_it,1) = dof(node,2);
                bc_it = bc_it + 1;
                bc(bc_it,1) = dof(node,3);
                bc_it = bc_it + 1;
            end
            
            % --- Fifth node ---
            node = node+z_step-y_step;
            edof(it,14:16) = dof(node,:);
            if el == 0
                bc(bc_it,1) = dof(node,1);
                bc_it = bc_it + 1;
                bc(bc_it,1) = dof(node,2);
                bc_it = bc_it + 1;
                bc(bc_it,1) = dof(node,3);
                bc_it = bc_it + 1;
            end
            
            % --- Sixth node ---
            node = node+x_step;
            edof(it,17:19) = dof(node,:);
            if el == nel_x-1
                bc(bc_it,1) = dof(node,1);
                bc_it = bc_it + 1;
                bc(bc_it,1) = dof(node,2);
                bc_it = bc_it + 1;
                bc(bc_it,1) = dof(node,3);
                bc_it = bc_it + 1;
            end
            
            % --- Seventh node ---
            node = node+y_step;
            edof(it,20:22) = dof(node,:);
            if el == nel_x-1
                bc(bc_it,1) = dof(node,1);
                bc_it = bc_it + 1;
                bc(bc_it,1) = dof(node,2);
                bc_it = bc_it + 1;
                bc(bc_it,1) = dof(node,3);
                bc_it = bc_it + 1;
            end
            
            % --- Eighth node ---
            node = node-x_step;
            edof(it,23:25) = dof(node,:);
            if el == 0
                bc(bc_it,1) = dof(node,1);
                bc_it = bc_it + 1;
                bc(bc_it,1) = dof(node,2);
                bc_it = bc_it + 1;
                bc(bc_it,1) = dof(node,3);
                bc_it = bc_it + 1;
            end

            
            
            % Reset node
            if el == nel_x-1 % If last element
                node = node-z_step-y_step+2;
            else % Otherwise, first node for next el. = second node for current
                node = node+x_step-y_step-z_step;
            end
            
            it = it+1;

        end
    end
end

% --- Creating global Stiffness & Force matrices ---

[ex,ey,ez] = coordxtr(edof,coord,dof,8);

ep = [2]; % No. integration points

% Material parameters for steel
E = 210000000; % Modulus of elasticity [Pa]
v = 0.3; % Poisson's ratio
D = hooke(4,E,v); % Material matrix

g=9.82; % Gravitational constant
rho = 7850; % Density [kg/m^3]
eq = [0; -g*rho; 0]; % Distibuted load vector [N/m^3]

f = zeros(ndof,1);
K = zeros(ndof);
for i=(1:nel) % Assembling
    [Ke,fe] = soli8e(ex(i,:), ey(i,:), ez(i,:), ep, D, eq);
    [K,f] = assem(edof(i,:), K, Ke, f, fe);
end

% --- Applying point forces  ---

point_force = -5000; % N
f(force_dofs) = f(force_dofs) + point_force;


% First a regular linear analysis is done, self-weight only

% --- Solving system of equations ---

a = solveq(K, f, bc);

% --- Extracting global displacements ---

ed = extract(edof,a);

% --- Extracting global displacements & calculating element stresses ---

es = zeros(ep*ep*ep,6,nel);
et = zeros(ep*ep*ep,6,nel);
eci = zeros(ep*ep*ep,3,nel);
for i=(1:nel)
    [es(:,:,i),et(:,:,i),eci(:,:,i)] = soli8s(ex(i,:),ey(i,:),ez(i,:),ep,D,ed(i,:));
end

% --- Calculating Principal stresses ---

%ps = zeros(3,nel);
%for i=(1:nel)
%    sigma_11 = mean(es(:,1,i));
%    sigma_22 = mean(es(:,2,i));
%    sigma_33 = mean(es(:,3,i));
%    sigma_12 = mean(es(:,4,i));
%    sigma_13 = mean(es(:,5,i));
%    sigma_23 = mean(es(:,6,i));

%    ps(1,i) = sigma_11 + sigma_22 + sigma_33;
%    ps(2,i) = sigma_11*sigma_22 + sigma_22*sigma_33 + sigma_11*sigma_33 - (sigma_12^2 + sigma_13^2 + sigma_23^2);
%    ps(3,i) = sigma_11*sigma_22*sigma_33 + 2*(sigma_12+sigma_13+sigma_23) - (sigma_12^2)*sigma_33 - (sigma_13^2)*sigma_22 - (sigma_23^2)*sigma_11;
%end

% --- Calculating Stress Tensors ---

Stress_tensors = zeros(3,3,nel);
for i=(1:nel)
    sigma_11 = mean(es(:,1,i));
    sigma_22 = mean(es(:,2,i));
    sigma_33 = mean(es(:,3,i));
    sigma_12 = mean(es(:,4,i));
    sigma_13 = mean(es(:,5,i));
    sigma_23 = mean(es(:,6,i));

    Stress_tensors(:,:,i) = [sigma_11 sigma_12 sigma_13; sigma_12 sigma_22 sigma_23; sigma_13 sigma_23 sigma_33];
end

% --- Calculating Principal stresses ---



% --- Using element stresses to calculate von Mises ---

vM_el = zeros(nel,1);
for i=(1:nel)
	%vM_el(i) = sqrt( 0.5 * ( (mean(es(:,1,i))-mean(es(:,2,i)))^2 + (mean(es(:,2,i))-mean(es(:,3,i)))^2 + (mean(es(:,3,i))-mean(es(:,1,i))) )^2 + 3 * ((mean(es(:,4,i)))^2 + (mean(es(:,5,i)))^2 + (mean(es(:,6,i)))^2) );
    s_xx = mean(es(:,1,i));
	s_yy = mean(es(:,2,i));
	s_zz = mean(es(:,3,i));
	s_xy = mean(es(:,4,i));
	s_xz = mean(es(:,5,i));
	s_yz = mean(es(:,6,i));
    vM_el(i) = sqrt( 0.5*((s_xx-s_yy)^2 + (s_yy-s_zz)^2 + (s_xx-s_yy)^2) + 3*(s_xy^2 + s_xz^2 + s_yz^2) );
end

% --- Nodal stresses for elements ---

%calc = (transpose(N) * N) \ transpose(N); % saving for quicker calculation

%ns = zeros(8,6,nel);
%for i = (1:nel)
%    ns(:,1,i) = calc*es(:,1,i);
%    ns(:,2,i) = calc*es(:,2,i);
%    ns(:,3,i) = calc*es(:,3,i);
%    ns(:,4,i) = calc*es(:,4,i);
%    ns(:,5,i) = calc*es(:,5,i);
%    ns(:,6,i) = calc*es(:,6,i);
%end

% --- Using element stresses to calculate von Mises ---

%vM_n = zeros(nel,8);
%for i=(1:nel)
%    s_xx = ns(:,1,i);
%	s_yy = ns(:,1,i);
%	s_zz = ns(:,1,i);
%	s_xy = ns(:,1,i);
%	s_xz = ns(:,1,i);
%	s_yz = ns(:,1,i);

%    vM_n(i,1) = sqrt(0.5 * ( (s_xx(1)-s_yy(1)).^2 + (s_yy(1)-s_zz(1)).^2 + (s_zz(1)-s_xx(1)).^2 ) + 3 * ((s_xy(1)).^2 + (s_xz(1)).^2 + (s_yz(1)).^2));
%    vM_n(i,2) = sqrt(0.5 * ( (s_xx(2)-s_yy(2)).^2 + (s_yy(2)-s_zz(2)).^2 + (s_zz(2)-s_xx(2)).^2 ) + 3 * ((s_xy(2)).^2 + (s_xz(2)).^2 + (s_yz(2)).^2));
%    vM_n(i,3) = sqrt(0.5 * ( (s_xx(3)-s_yy(3)).^2 + (s_yy(3)-s_zz(3)).^2 + (s_zz(3)-s_xx(3)).^2 ) + 3 * ((s_xy(3)).^2 + (s_xz(3)).^2 + (s_yz(3)).^2));
%    vM_n(i,4) = sqrt(0.5 * ( (s_xx(4)-s_yy(4)).^2 + (s_yy(4)-s_zz(4)).^2 + (s_zz(4)-s_xx(4)).^2 ) + 3 * ((s_xy(4)).^2 + (s_xz(4)).^2 + (s_yz(4)).^2));
%    vM_n(i,5) = sqrt(0.5 * ( (s_xx(5)-s_yy(5)).^2 + (s_yy(5)-s_zz(5)).^2 + (s_zz(5)-s_xx(5)).^2 ) + 3 * ((s_xy(5)).^2 + (s_xz(5)).^2 + (s_yz(5)).^2));
%    vM_n(i,6) = sqrt(0.5 * ( (s_xx(6)-s_yy(6)).^2 + (s_yy(6)-s_zz(6)).^2 + (s_zz(6)-s_xx(6)).^2 ) + 3 * ((s_xy(6)).^2 + (s_xz(6)).^2 + (s_yz(6)).^2));
%    vM_n(i,7) = sqrt(0.5 * ( (s_xx(7)-s_yy(7)).^2 + (s_yy(7)-s_zz(7)).^2 + (s_zz(7)-s_xx(7)).^2 ) + 3 * ((s_xy(7)).^2 + (s_xz(7)).^2 + (s_yz(7)).^2));
%    vM_n(i,8) = sqrt(0.5 * ( (s_xx(8)-s_yy(8)).^2 + (s_yy(8)-s_zz(8)).^2 + (s_zz(8)-s_xx(8)).^2 ) + 3 * ((s_xy(8)).^2 + (s_xz(8)).^2 + (s_yz(8)).^2));

    %vM_n(i,1) = sqrt( 0.5 * ( (ns(1,1,i)-ns(1,2,i)).^2 + (ns(1,2,i)-ns(1,3,i)).^2 + (ns(1,3,i)-ns(1,1,i)).^2 ) + 3 * ((ns(1,4,i)).^2 + (ns(1,5,i)).^2 + (ns(1,6,i)).^2) );
    %vM_n(i,2) = sqrt( 0.5 * ( (ns(2,1,i)-ns(2,2,i)).^2 + (ns(2,2,i)-ns(2,3,i)).^2 + (ns(2,3,i)-ns(2,1,i)).^2 ) + 3 * ((ns(2,4,i)).^2 + (ns(2,5,i)).^2 + (ns(2,6,i)).^2) );
    %vM_n(i,3) = sqrt( 0.5 * ( (ns(3,1,i)-ns(3,2,i)).^2 + (ns(3,2,i)-ns(3,3,i)).^2 + (ns(3,3,i)-ns(3,1,i)).^2 ) + 3 * ((ns(3,4,i)).^2 + (ns(3,5,i)).^2 + (ns(3,6,i)).^2) );
    %vM_n(i,4) = sqrt( 0.5 * ( (ns(4,1,i)-ns(4,2,i)).^2 + (ns(4,2,i)-ns(4,3,i)).^2 + (ns(4,3,i)-ns(4,1,i)).^2 ) + 3 * ((ns(4,4,i)).^2 + (ns(4,5,i)).^2 + (ns(4,6,i)).^2) );
    %vM_n(i,5) = sqrt( 0.5 * ( (ns(5,1,i)-ns(5,2,i)).^2 + (ns(5,2,i)-ns(5,3,i)).^2 + (ns(5,3,i)-ns(5,1,i)).^2 ) + 3 * ((ns(5,4,i)).^2 + (ns(5,5,i)).^2 + (ns(5,6,i)).^2) );
    %vM_n(i,6) = sqrt( 0.5 * ( (ns(6,1,i)-ns(6,2,i)).^2 + (ns(6,2,i)-ns(6,3,i)).^2 + (ns(6,3,i)-ns(6,1,i)).^2 ) + 3 * ((ns(6,4,i)).^2 + (ns(6,5,i)).^2 + (ns(6,6,i)).^2) );
    %vM_n(i,7) = sqrt( 0.5 * ( (ns(7,1,i)-ns(7,2,i)).^2 + (ns(7,2,i)-ns(7,3,i)).^2 + (ns(7,3,i)-ns(7,1,i)).^2 ) + 3 * ((ns(7,4,i)).^2 + (ns(7,5,i)).^2 + (ns(7,6,i)).^2) );
    %vM_n(i,8) = sqrt( 0.5 * ( (ns(8,1,i)-ns(8,2,i)).^2 + (ns(8,2,i)-ns(8,3,i)).^2 + (ns(8,3,i)-ns(8,1,i)).^2 ) + 3 * ((ns(8,4,i)).^2 + (ns(8,5,i)).^2 + (ns(8,6,i)).^2) );
%end

% % The nodal stresses are interpolated through and averaged in order to
% % avoid discontinuities in the visualization of the model
% 
% % Elements 1-50: bottom left
% % Elements 51-100: above elements 1-50
% % Elements 101-150: above elements 51-100
% % Elements 151-200: above elements 101-150
% 
% % Elements 201-250: bottom middle
% % Elements 251-300: above elements 201-250
% % Elements 301-350: above elements 251-300
% % Elements 351-400: above elements 301-350
% 
% % Elements 401-450: bottom right
% % Elements 451-500: above elements 401-450
% % Elements 501-550: above elements 451-500
% % Elements 551-600: above elements 501-550
% 
% x = 1; % Next element in x-direction
% y = nel_x; % Next element in y-direction
% z = y*nel_y; % Next element in z-direction
% 
% %x_step = 1; % Next node in x-direction
% %y_step = nnode_x; % Next node in y-direction
% %z_step = (y_step)*nnode_y; % Next node in z-direction
% 
% el = zeros(3,3,3); % Used to store adjacent elements, 26+1 where (2,2,2) is the studied element
% 
% ns = zeros(8,6,nel);
% for i = (1:nel)
%     %elements(2,2,2) = i;
% 
%     el(1,1,1) = i -x -y -z;
%     el(2,1,1) = i    -y -z;
%     el(3,1,1) = i +x -y -z;
% 
%     el(1,2,1) = i -x -z;
%     el(2,2,1) = i    -z;
%     el(3,2,1) = i +x -z;
% 
%     el(1,3,1) = i -x +y -z;
%     el(2,3,1) = i    +y -z;
%     el(3,3,1) = i +x +y -z;
% 
%     el(1,1,2) = i -x -y;
%     el(2,1,2) = i    -y;
%     el(3,1,2) = i +x -y;
% 
%     el(1,2,2) = i -x;
%     el(2,2,2) = i;
%     el(3,2,2) = i +x;
% 
%     el(1,3,2) = i -x +y -z;
%     el(2,3,2) = i    +y-z;
%     el(3,3,2) = i +x +y -z;
% 
%     el(1,1,3) = i -x -y +z;
%     el(2,1,3) = i    -y +z;
%     el(3,1,3) = i +x -y +z;
% 
%     el(1,2,3) = i -x +z;
%     el(2,2,3) = i    +z;
%     el(3,2,3) = i +x +z;
% 
%     el(1,3,3) = i -x +y +z;
%     el(2,3,3) = i    +y +z;
%     el(3,3,3) = i +x +y +z;
% 
%     el(el<0) = 0; % if an element is negative (doesn't exist), set to 0
%     
%     %if el(1,1,1) == 0 %
%     %    el(1,1,2) = 0;
%     %end
% 
%     if el(1,2,1) == 0 %
%         el(1,1,1) = 0;
%         el(1,3,1) = 0;
% 
%         el(1,1,2) = 0;
%         el(1,2,2) = 0;
%         el(1,3,2) = 0;
% 
%         el(1,1,3) = 0;
%         el(1,2,3) = 0;
%         el(1,3,3) = 0;
%     end
% 
%     if el(2,1,1) == 0 %
%         el(1,1,1) = 0;
%         el(3,1,1) = 0;
% 
%         el(1,1,2) = 0;
%         el(2,1,2) = 0;
%         el(3,1,2) = 0;
% 
%         el(1,1,3) = 0;
%         el(2,1,3) = 0;
%         el(3,1,3) = 0;
%     end
% 
%     %if el(1,3,1) == 0 %
%     %    el(1,3,1) = 0;
%     %end
% 
%     % Corner -> No interpolation
%     % Only elements in x-direction -> Interpolated with one value in x-dir.
%     % Only elements in y-direction -> Interpolated with one value in y-dir.
%     % Only elements in z-direction -> Interpolated with one value in z-dir.
%     % No element behind/ahead  -> Interpolated with value in y- & z-dir.
%     % No element below/above  -> Interpolated with value in x- & z-dir.
%     % No element to left/right  -> Interpolated with value in x- & y-dir.
% 
%     % Behind/ahead refers to -/+ x-direction
%     % Below/above refers to -/+ y-direction
%     % Left/right refers to -/+ z-direction
% 
%     for j = (1:8)
%         if j==1 % Node 1 for element
%             if el(1,2,2) == 0 && el(2,1,2) == 0 && el(2,2,1) == 0 % If corner
%                 %disp('corner')
%                 ns(j,:,i) = ns_old(1,:,el(2,2,2));
%             elseif el(1,2,2) == 0 && el(2,1,2) == 0 % If only elements in z-direction
%                 %disp('only z')
%                 ns(j,:,i) = mean( [ns_old(1,:,el(2,2,2)), ...
%                     ns_old(4,:,el(2,2,1))] );
%             elseif el(1,2,2) == 0 && el(2,2,1) == 0 % If only elements in y-direction
%                 %disp('only y')
%                 ns(j,:,i) = mean( [ns_old(1,:,el(2,2,2)), ...
%                     ns_old(5,:,el(2,1,2))] );
%             elseif el(2,1,2) == 0 && el(2,2,1) == 0 % If only elements in x-direction
%                 %disp('only x')
%                 ns(j,:,i) = mean( [ns_old(1,:,el(2,2,2)), ...
%                     ns_old(2,:,el(1,2,2))] );
%             elseif el(1,2,2) == 0 % If no elements behind (x-direction)
%                 %disp('no x')
%                 ns(j,:,i) = mean( [ns_old(1,:,el(2,2,2)), ...
%                     ns_old(4,:,el(2,2,1)), ...
%                     ns_old(5,:,el(2,1,2)), ...
%                     ns_old(8,:,el(2,1,1))] );
%             elseif el(2,1,2) == 0 % If no elements below (y-direction)
%                 %disp('no y')
%                 ns(j,:,i) = mean( [ns_old(1,:,el(2,2,2)), ...
%                     ns_old(2,:,el(1,2,2)), ...
%                     ns_old(4,:,el(2,2,1)), ...
%                     ns_old(6,:,el(1,2,1))] );
%             elseif el(2,2,1) == 0 % If no elements to the left (z-direction)
%                 %disp('no z')
%                 ns(j,:,i) = mean( [ns_old(1,:,el(2,2,2)), ...
%                     ns_old(2,:,el(1,2,2)), ...
%                     ns_old(3,:,el(1,1,2)), ...
%                     ns_old(5,:,el(2,1,2))] );
%             else % Otherwise, elements in all directions
%                 %disp('all dir')
%                 ns(j,:,i) = mean( [ns_old(1,:,el(2,2,2)), ...
%                     ns_old(2,:,el(1,2,2)), ...
%                     ns_old(3,:,el(1,1,2)), ...
%                     ns_old(4,:,el(2,2,1)), ...
%                     ns_old(5,:,el(2,1,2)), ...
%                     ns_old(6,:,el(1,2,1)), ...
%                     ns_old(7,:,el(1,1,1)), ...
%                     ns_old(8,:,el(2,1,1))] );
%             end
%         
% 
% 
% 
% 
% 
%         elseif j==2 % Node 2 for element
%             if el(3,2,2) == 0 && el(2,1,2) == 0 && el(2,2,1) == 0 % If corner
%                 %disp('corner')
%                 ns(j,:,i) = ns_old(2,:,el(2,2,2));
%             elseif el(3,2,2) == 0 && el(2,1,2) == 0 % If only elements in z-direction
%                 %disp('only z')
%                 ns(j,:,i) = mean( [ns_old(2,:,el(2,2,2)), ...
%                     ns_old(6,:,el(2,2,1))] );
%             elseif el(3,2,2) == 0 && el(2,2,1) == 0 % If only elements in y-direction
%                 %disp('only y')
%                 ns(j,:,i) = mean( [ns_old(2,:,el(2,2,2)), ...
%                     ns_old(3,:,el(2,1,2))]  );
%             elseif el(2,1,2) == 0 && el(2,2,1) == 0 % If only elements in x-direction
%                 %disp('only x')
%                 ns(j,:,i) = mean( [ns_old(2,:,el(2,2,2)), ...
%                     ns_old(1,:,el(3,2,2))] );
%             elseif el(3,2,2) == 0 % If no elements ahead (x-direction)
%                 %disp('no x')
%                 ns(j,:,i) = mean( [ns_old(2,:,el(2,2,2)), ...
%                     ns_old(6,:,el(2,2,1)), ...
%                     ns_old(3,:,el(2,1,2)), ...
%                     ns_old(7,:,el(2,1,1))] );
%             elseif el(2,1,2) == 0 % If no elements below (y-direction)
%                 %disp('no y')
%                 ns(j,:,i) = mean( [ns_old(2,:,el(2,2,2)), ...
%                     ns_old(1,:,el(3,2,2)), ...
%                     ns_old(6,:,el(2,2,1)), ...
%                     ns_old(5,:,el(3,2,1))] );
%             elseif el(2,2,1) == 0 % If no elements to the left (z-direction)
%                 %disp('no z')
%                 ns(j,:,i) = mean( [ns_old(2,:,el(2,2,2)), ...
%                     ns_old(1,:,el(3,2,2)), ...
%                     ns_old(4,:,el(3,1,2)), ...
%                     ns_old(3,:,el(2,1,2))] );
%             else % Otherwise, elements in all directions
%                 %disp('all dir')
%                 ns(j,:,i) = mean( [ns_old(2,:,el(2,2,2)), ...
%                     ns_old(1,:,el(3,2,2)), ...
%                     ns_old(4,:,el(3,1,2)), ...
%                     ns_old(6,:,el(2,2,1)), ...
%                     ns_old(3,:,el(2,1,2)), ...
%                     ns_old(5,:,el(3,2,1)), ...
%                     ns_old(8,:,el(3,1,1)), ...
%                     ns_old(7,:,el(2,1,1))] );
%             end
%         
% 
% 
%         
% 
% 
% 
%         elseif j==3 % Node 3 for element
%             if el(3,2,2) == 0 && el(2,3,2) == 0 && el(2,2,1) == 0 % If corner
%                 %disp('corner')
%                 ns(j,:,i) = ns_old(3,:,el(2,2,2));
%             elseif el(3,2,2) == 0 && el(2,3,2) == 0 % If only elements in z-direction
%                 %disp('only z')
%                 ns(j,:,i) = mean( [ns_old(3,:,el(2,2,2)), ...
%                     ns_old(6,:,el(2,2,1))] );
%             elseif el(3,2,2) == 0 && el(2,2,1) == 0 % If only elements in y-direction
%                 %disp('only y')
%                 ns(j,:,i) = mean( ns_old(3,:,el(2,2,2)), ...
%                     + ns_old(2,:,el(2,3,2)) );
%             elseif el(2,3,2) == 0 && el(2,2,1) == 0 % If only elements in x-direction
%                 %disp('only x')
%                 ns(j,:,i) = mean( ns_old(3,:,el(2,2,2)) ...
%                     + ns_old(4,:,el(3,2,2)) );
%             elseif el(3,2,2) == 0 % If no elements ahead (x-direction)
%                 %disp('no x')
%                 ns(j,:,i) = mean( ns_old(3,:,el(2,2,2)) ...
%                     + ns_old(6,:,el(2,2,1)) ...
%                     + ns_old(2,:,el(2,3,2)) ...
%                     + ns_old(7,:,el(2,3,1)) );
%             elseif el(2,3,2) == 0 % If no elements above (y-direction)
%                 %disp('no y')
%                 ns(j,:,i) = mean( ns_old(3,:,el(2,2,2)) ...
%                     + ns_old(4,:,el(3,2,2)) ...
%                     + ns_old(6,:,el(2,2,1)) ...
%                     + ns_old(8,:,el(3,2,1)) );
%             elseif el(2,2,1) == 0 % If no elements to the left (z-direction)
%                 ns(j,:,i) = mean( ns_old(3,:,el(2,2,2)) ...
%                     + ns_old(4,:,el(3,2,2)) ...
%                     + ns_old(1,:,el(3,3,2)) ...
%                     + ns_old(2,:,el(2,3,2)) );
%             else % Otherwise, elements in all directions
%                 %disp('all dir')
%                 ns(j,:,i) = mean( ns_old(3,:,el(2,2,2)) ...
%                     + ns_old(4,:,el(3,2,2)) ...
%                     + ns_old(1,:,el(3,3,2)) ...
%                     + ns_old(6,:,el(2,2,1)) ...
%                     + ns_old(2,:,el(2,3,2)) ...
%                     + ns_old(8,:,el(3,2,1)) ...
%                     + ns_old(5,:,el(3,3,1)) ...
%                     + ns_old(7,:,el(2,3,1)) );
%             end
%         
% 
% 
% 
%         
% 
% 
% 
% 
% 
%         elseif j==4 % Node 4 for element
%             if el(1,2,2) == 0 && el(2,3,2) == 0 && el(2,2,1) == 0 % If corner
%                 %disp('corner')
%                 ns(j,:,i) = ns_old(4,:,el(2,2,2));
%             elseif el(1,2,2) == 0 && el(2,3,2) == 0 % If only elements in z-direction
%                 %disp('only z')
%                 ns(j,:,i) = mean( ns_old(4,:,el(2,2,2)) ...
%                     + ns_old(8,:,el(2,2,1)) );
%             elseif el(1,2,2) == 0 && el(2,2,1) == 0 % If only elements in y-direction
%                 %disp('only y')
%                 ns(j,:,i) = mean( ns_old(4,:,el(2,2,2)) ...
%                     + ns_old(1,:,el(2,3,2)) );
%             elseif el(2,3,2) == 0 && el(2,2,1) == 0 % If only elements in x-direction
%                 %disp('only x')
%                 ns(j,:,i) = mean( ns_old(4,:,el(2,2,2)) ...
%                     + ns_old(3,:,el(1,2,2)) );
%             elseif el(1,2,2) == 0 % If no elements behind (x-direction)
%                 %disp('no x')
%                 ns(j,:,i) = mean( ns_old(4,:,el(2,2,2)) ...
%                     + ns_old(8,:,el(2,2,1)) ...
%                     + ns_old(1,:,el(2,3,2)) ...
%                     + ns_old(5,:,el(2,3,1)) );
%             elseif el(2,3,2) == 0 % If no elements above (y-direction)
%                 %disp('no y')
%                 ns(j,:,i) = mean( ns_old(4,:,el(2,2,2)) ...
%                     + ns_old(3,:,el(1,2,2)) ...
%                     + ns_old(8,:,el(2,2,1)) ...
%                     + ns_old(7,:,el(1,2,1)) );
%             elseif el(2,2,1) == 0 % If no elements to the left (z-direction)
%                 ns(j,:,i) = mean( ns_old(4,:,el(2,2,2)) ...
%                     + ns_old(3,:,el(1,2,2)) ...
%                     + ns_old(2,:,el(1,3,2)) ...
%                     + ns_old(1,:,el(2,3,2)) );
%             else % Otherwise, elements in all directions
%                 %disp('all dir')
%                 ns(j,:,i) = mean( ns_old(4,:,el(2,2,2)) ...
%                     + ns_old(3,:,el(1,2,2)) ...
%                     + ns_old(2,:,el(1,3,2)) ...
%                     + ns_old(8,:,el(2,2,1)) ...
%                     + ns_old(1,:,el(2,3,2)) ...
%                     + ns_old(7,:,el(1,2,1)) ...
%                     + ns_old(6,:,el(1,3,1)) ...
%                     + ns_old(5,:,el(2,3,1)) );
%             end
% 
% 
% 
% 
% 
%         elseif j==5 % Node 5 for element
%             if el(1,2,2) == 0 && el(2,1,2) == 0 && el(2,2,3) == 0 % If corner
%                 %disp('corner')
%                 ns(j,:,i) = ns_old(5,:,el(2,2,2));
%             elseif el(1,2,2) == 0 && el(2,1,2) == 0 % If only elements in z-direction
%                 %disp('only z')
%                 ns(j,:,i) = mean( ns_old(5,:,el(2,2,2)) ...
%                     + ns_old(1,:,el(2,2,3)) );
%             elseif el(1,2,2) == 0 && el(2,2,3) == 0 % If only elements in y-direction
%                 %disp('only y')
%                 ns(j,:,i) = mean( ns_old(5,:,el(2,2,2)) ...
%                     + ns_old(8,:,el(2,1,2)) );
%             elseif el(2,1,2) == 0 && el(2,2,3) == 0 % If only elements in x-direction
%                 %disp('only x')
%                 ns(j,:,i) = mean( ns_old(5,:,el(2,2,2)) ...
%                     + ns_old(6,:,el(1,2,2)) );
%             elseif el(1,2,2) == 0 % If no elements behind (x-direction)
%                 %disp('no x')
%                 ns(j,:,i) = mean( ns_old(5,:,el(2,2,2)) ...
%                     + ns_old(1,:,el(2,2,3)) ...
%                     + ns_old(8,:,el(2,1,2)) ...
%                     + ns_old(4,:,el(2,1,3)) );
%             elseif el(2,1,2) == 0 % If no elements below (y-direction)
%                 %disp('no y')
%                 ns(j,:,i) = mean( ns_old(5,:,el(2,2,2)) ...
%                     + ns_old(6,:,el(1,2,2)) ...
%                     + ns_old(1,:,el(2,2,3)) ...
%                     + ns_old(2,:,el(1,2,3)) );
%             elseif el(2,2,3) == 0 % If no elements to the right (z-direction)
%                 %disp('no z')
%                 ns(j,:,i) = mean( ns_old(5,:,el(2,2,2)) ...
%                     + ns_old(6,:,el(1,2,2)) ...
%                     + ns_old(7,:,el(1,1,2)) ...
%                     + ns_old(8,:,el(2,1,2)) );
%             else % Otherwise, elements in all directions
%                 %disp('all dir')
%                 ns(j,:,i) = mean( ns_old(5,:,el(2,2,2)) ...
%                     + ns_old(6,:,el(1,2,2)) ...
%                     + ns_old(7,:,el(1,1,2)) ...
%                     + ns_old(1,:,el(2,2,3)) ...
%                     + ns_old(8,:,el(2,1,2)) ...
%                     + ns_old(2,:,el(1,2,3)) ...
%                     + ns_old(3,:,el(1,1,3)) ...
%                     + ns_old(4,:,el(2,1,3)) );
%             end
%         
%         
%         elseif j==6 % Node 6 for element
%             if el(3,2,2) == 0 && el(2,1,2) == 0 && el(2,2,3) == 0 % If corner
%                 %disp('corner')
%                 ns(j,:,i) = ns_old(6,:,el(2,2,2));
%             elseif el(3,2,2) == 0 && el(2,1,2) == 0 % If only elements in z-direction
%                 %disp('only z')
%                 ns(j,:,i) = mean( ns_old(6,:,el(2,2,2)) ...
%                     + ns_old(2,:,el(2,2,3)) );
%             elseif el(3,2,2) == 0 && el(2,2,3) == 0 % If only elements in y-direction
%                 %disp('only y')
%                 ns(j,:,i) = mean( ns_old(6,:,el(2,2,2)) ...
%                     + ns_old(7,:,el(2,1,2)) );
%             elseif el(2,1,2) == 0 && el(2,2,3) == 0 % If only elements in x-direction
%                 %disp('only x')
%                 ns(j,:,i) = mean( ns_old(6,:,el(2,2,2)) ...
%                     + ns_old(5,:,el(3,2,2)) );
%             elseif el(3,2,2) == 0 % If no elements ahead (x-direction)
%                 %disp('no x')
%                 ns(j,:,i) = mean( ns_old(6,:,el(2,2,2)) ...
%                     + ns_old(2,:,el(2,2,3)) ...
%                     + ns_old(7,:,el(2,1,2)) ...
%                     + ns_old(3,:,el(2,1,3)) );
%             elseif el(2,1,2) == 0 % If no elements below (y-direction)
%                 %disp('no y')
%                 ns(j,:,i) = mean( ns_old(6,:,el(2,2,2)) ...
%                     + ns_old(5,:,el(3,2,2)) ...
%                     + ns_old(2,:,el(2,2,3)) ...
%                     + ns_old(1,:,el(3,2,3)) );
%             elseif el(2,2,3) == 0 % If no elements to the right (z-direction)
%                 %disp('no z')
%                 ns(j,:,i) = mean( ns_old(6,:,el(2,2,2)) ...
%                     + ns_old(5,:,el(3,2,2)) ...
%                     + ns_old(8,:,el(3,1,2)) ...
%                     + ns_old(7,:,el(2,1,2)) );
%             else % Otherwise, elements in all directions
%                 %disp('all dir')
%                 ns(j,:,i) = mean( ns_old(6,:,el(2,2,2)) ...
%                     + ns_old(5,:,el(3,2,2)) ...
%                     + ns_old(8,:,el(3,1,2)) ...
%                     + ns_old(2,:,el(2,2,3)) ...
%                     + ns_old(7,:,el(2,1,2)) ...
%                     + ns_old(1,:,el(3,2,3)) ...
%                     + ns_old(4,:,el(3,1,3)) ...
%                     + ns_old(3,:,el(2,1,3)) );
%             end
%         
%         
%         elseif j==7 % Node 7 for element
%             if el(3,2,2) == 0 && el(2,3,2) == 0 && el(2,2,3) == 0 % If corner
%                 %disp('corner')
%                 ns(j,:,i) = ns_old(7,:,el(2,2,2));
%             elseif el(3,2,2) == 0 && el(2,3,2) == 0 % If only elements in z-direction
%                 %disp('only z')
%                 ns(j,:,i) = mean( ns_old(7,:,el(2,2,2)) ...
%                     + ns_old(3,:,el(2,2,3)) );
%             elseif el(3,2,2) == 0 && el(2,2,3) == 0 % If only elements in y-direction
%                 %disp('only y')
%                 ns(j,:,i) = mean( ns_old(7,:,el(2,2,2)) ...
%                     + ns_old(6,:,el(2,3,2)) );
%             elseif el(2,3,2) == 0 && el(2,2,3) == 0 % If only elements in x-direction
%                 %disp('only x')
%                 ns(j,:,i) = mean( ns_old(7,:,el(2,2,2)) ...
%                     + ns_old(8,:,el(3,2,2)) );
%             elseif el(3,2,2) == 0 % If no elements ahead (x-direction)
%                 %disp('no x')
%                 ns(j,:,i) = mean( ns_old(7,:,el(2,2,2)) ...
%                     + ns_old(3,:,el(2,2,3)) ...
%                     + ns_old(6,:,el(2,3,2)) ...
%                     + ns_old(2,:,el(2,3,3)) );
%             elseif el(2,3,2) == 0 % If no elements above (y-direction)
%                 %disp('no y')
%                 ns(j,:,i) = mean( ns_old(7,:,el(2,2,2)) ...
%                     + ns_old(8,:,el(3,2,2)) ...
%                     + ns_old(3,:,el(2,2,3)) ...
%                     + ns_old(4,:,el(3,2,3)) );
%             elseif el(2,2,3) == 0 % If no elements to the right (z-direction)
%                 %disp('no z')
%                 ns(j,:,i) = mean( ns_old(7,:,el(2,2,2)) ...
%                     + ns_old(8,:,el(3,2,2)) ...
%                     + ns_old(5,:,el(3,3,2)) ...
%                     + ns_old(6,:,el(2,3,2)) );
%             else % Otherwise, elements in all directions
%                 %disp('all dir')
%                 ns(j,:,i) = mean( ns_old(7,:,el(2,2,2)) ...
%                     + ns_old(8,:,el(3,2,2)) ...
%                     + ns_old(5,:,el(3,3,2)) ...
%                     + ns_old(3,:,el(2,2,3)) ...
%                     + ns_old(6,:,el(2,3,2)) ...
%                     + ns_old(4,:,el(3,2,3)) ...
%                     + ns_old(1,:,el(3,3,3)) ...
%                     + ns_old(2,:,el(2,3,3)) );
%             end
%         
%         
%         elseif j==8 % Node 8 for element
%             if el(1,2,2) == 0 && el(2,3,2) == 0 && el(2,2,3) == 0 % If corner
%                 %disp('corner')
%                 ns(j,:,i) = ns_old(8,:,el(2,2,2));
%             elseif el(1,2,2) == 0 && el(2,3,2) == 0 % If only elements in z-direction
%                 %disp('only z')
%                 ns(j,:,i) = mean( [ns_old(8,:,el(2,2,2)), ...
%                     ns_old(4,:,el(2,2,3))] );
%             elseif el(1,2,2) == 0 && el(2,2,3) == 0 % If only elements in y-direction
%                 %disp('only y')
%                 ns(j,:,i) = mean( [ns_old(8,:,el(2,2,2)), ...
%                     ns_old(5,:,el(2,3,2))] );
%             elseif el(2,3,2) == 0 && el(2,2,3) == 0 % If only elements in x-direction
%                 %disp('only x')
%                 ns(j,:,i) = mean( [ns_old(8,:,el(2,2,2)), ...
%                     ns_old(7,:,el(1,2,2))] );
%             elseif el(1,2,2) == 0 % If no elements behind (x-direction)
%                 %disp('no x')
%                 ns(j,:,i) = mean( [ns_old(8,:,el(2,2,2)), ...
%                     ns_old(4,:,el(2,2,3)), ...
%                     ns_old(5,:,el(2,3,2)), ...
%                     ns_old(1,:,el(2,3,3))] );
%             elseif el(2,3,2) == 0 % If no elements above (y-direction)
%                 %disp('no y')
%                 ns(j,:,i) = mean( [ns_old(8,:,el(2,2,2)), ...
%                     ns_old(7,:,el(1,2,2)), ...
%                     ns_old(4,:,el(2,2,3)), ...
%                     ns_old(3,:,el(1,2,3))] );
%             elseif el(2,2,3) == 0 % If no elements to the right (z-direction)
%                 %disp('no z')
%                 ns(j,:,i) = mean( [ns_old(8,:,el(2,2,2)), ...
%                     ns_old(7,:,el(1,2,2)), ...
%                     ns_old(6,:,el(1,3,2)), ...
%                     ns_old(5,:,el(2,3,2))] );
%             else % Otherwise, elements in all directions
%                 %disp('all dir')
%                 ns(j,:,i) = mean( [ns_old(8,:,el(2,2,2)) ...
%                     ns_old(7,:,el(1,2,2)), ...
%                     ns_old(6,:,el(1,3,2)), ...
%                     ns_old(4,:,el(2,2,3)), ...
%                     ns_old(5,:,el(2,3,2)), ...
%                     ns_old(3,:,el(1,2,3)), ...
%                     ns_old(2,:,el(1,3,3)), ...
%                     ns_old(1,:,el(2,3,3))] );
%             end
%         
%         
%         
%         end
%     
% 
%     end
%     el = zeros(3,3,3);
% 
% 
% 
%     %el_behind = i-xstep;
%     %el_behind_below = i-xstep
% 
% 
% 
% 
% 
%     %el_infront = i+xstep;
% 
% 
%     %el_below = i-ystep;
%     %el_above = i+ystep;
% 
%     %el_left = i-zstep;
%     %el_right = i+zstep;
% 
%     
% 
%     
% 
%     % Coord 1, -x, -y, -z
% 
%     %if el_behind > 0 && el_below > 0 && el_left > 0
%     %    ns(1,:,i)
%     %elseif 
%     %else
%     %    ns(1,:,i) = ns_old(1,:,i);
%     %end
% 
%     
%     el = zeros(3,3,3);
% end



% Now an eigenvalue analysis of the model is done

% --- Masses for element ---

m = zeros(8);
for i = (1:8)
    for j = (1:8)
        m(i,j) = (rho*d*d*d/8)*(1+(1/3)*xsi(i)*xsi(j))*(1+(1/3)*eta(i)*eta(j))*(1+(1/3)*zet(i)*zet(j));
    end
end

% --- Element mass matrix ---

iter_i = 1;
iter_j = 1;

Me = zeros(3*8);
for i = (1:3:3*8)
    iter_j = 1;
    for j = (1:3:3*8)
        Me(i,j) = m(iter_i,iter_j);
        Me(i+1,j+1) = m(iter_i,iter_j);
        Me(i+2,j+2) = m(iter_i,iter_j);
        iter_j = iter_j+1;
    end
    iter_i = iter_i+1;
end

% --- Global mass matrix ---

M = zeros(ndof);
for i=(1:nel)
    M = assem(edof(i,:), M, Me);
end

% --- Eigenvalue analysis ---

[lambda,eig] = eigen(K,M,bc(:,1));



% --- Results from both analyses are saved in a .mat-file ---

save('exv4.mat','coord','dof','edof','bc','force_dofs','a','ed','Stress_tensors','vM_el','vM_n','lambda','eig')
