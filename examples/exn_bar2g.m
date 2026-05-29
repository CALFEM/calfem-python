% example exn_bar2g 
%--------------------------------------------------------------------------
% PURPOSE 
%    Analysis of a plane truss using second order theory.
%--------------------------------------------------------------------------

% LAST MODIFIED: O Dahlblom  2016-01-18
% Copyright (c)  Division of Structural Mechanics and
%                Division of Solid Mechanics.
%                Lund University
%--------------------------------------------------------------------------
 echo on

% ----- Topology -----

 Edof=[1  1  2  5  6;
       2  3  4  5  6]; 

% ----- Element properties and global coordinates ------------------------- 
      
 E=10e9;  
 A1=4e-2;     A2=1e-2;
 ep1=[E A1];  ep2=[E A2];

 ex1=[0 1.6]; ey1=[0 0];
 ex2=[0 1.6]; ey2=[1.2 0];

% ----- Initial values for the iteration ----------------------------------

 eps=1e-6;           % Error norm
 QX1=0.01;    QX2=0; % Initial axial forces
 QX01=1;             % Axial force of the initial former iteration
 n=0;                % Iteration counter

% ----- Iteration procedure -----------------------------------------------
 
 while(abs((QX1-QX01)/QX01)>eps)
 n=n+1

   K=zeros(6,6);
   f=zeros(6,1);	
   f(5)=-10e6;
   f(6)=-0.2e6;

   Ke1=bar2ge(ex1,ey1,ep1,QX1);
   Ke2=bar2ge(ex2,ey2,ep2,QX2);
   K=assem(Edof(1,:),K,Ke1);
   K=assem(Edof(2,:),K,Ke2);
   bc=[1 0;2 0;3 0;4 0];	
   [a,r]=solveq(K,f,bc)

   Ed=extract_ed(Edof,a);

   QX01=QX1; 
   [es1,QX1]=bar2gs(ex1,ey1,ep1,Ed(1,:))
   [es2,QX2]=bar2gs(ex2,ey2,ep2,Ed(2,:))
  
   if(n>20)
     disp('The solution does not converge')
     return
   end
 end
 echo off


