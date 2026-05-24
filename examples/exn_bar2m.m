% example exn_bar2m
%----------------------------------------------------------------
% PURPOSE 
%    Analysis of a plane truss with material nonlinearity.
%----------------------------------------------------------------

% LAST MODIFIED: O Dahlblom  2015-12-02
% Copyright (c)  Division of Structural Mechanics and
%                Division of Solid Mechanics.
%                Lund University
%----------------------------------------------------------------

echo off

edof=[1 1 2 5 6;
      2 5 6 7 8;
      3 3 4 5 6]

bc=[1 0; 2 0; 3 0; 4 0; 7 0; 8 0];

ex=[0 1.6; 1.6 1.6; 0 1.6]
ey=[0 0; 0 1.2; 1.2 0]

Em=200.0e9;
E=ones(3,1)*Em;
A1=6.0e-4; A2=3.0e-4; A3=10.0e-4;
A=[A1 A2 A3]'
SY=400.0e6;
Ns=SY*A

dp=4.0e3

incr=100;

a=zeros(8,1);
r=zeros(8,1);
es=zeros(3,1);

plbar=0;
pl(1,:)=[0 0];

% Forward Euler incremental solution

for i=1:incr
   K=zeros(8);
   df=zeros(8,1);
   df(6)=-dp;
  
   % Create and assemble element tangent stiffness matrices
   for j=1:3  
      ep=[E(j),A(j)];
      Ke=bar2e(ex(j,:),ey(j,:),ep);
      K=assem(edof(j,:),K,Ke);
   end;
      
   % Stop iteration if determinant det(Kr)<=0
   Kr=red(K,bc(:,1));
   if det(Kr)<=0 
      disp(['Determinant zero after increment ',num2str(i-1)])
      break;
  end;
   
   % Solve for the displacement increment and determine total displacements
   [da,dr]=solveq(K,df,bc);
   a=a+da;
   r=r+dr;
   
   % Determine normal forces in elements
   ded=extract_ed(edof,da);
   for j=1:3
     ep=[E(j),A(j)];
     desj=bar2s(ex(j,:),ey(j,:),ep,ded(j,:));
     des(j,1)=desj(1);
   end;
   es=es+des;
   for j=1:3
     E(j)=Em; if abs(es(j))>=Ns(j);  E(j)=0; end 
   end; 
   
   % Determine if the stress in a bar has reached the yield stress
   newplbar=sum(abs(es)>Ns);
   if newplbar > plbar
     plbar=newplbar;
     disp([num2str(plbar),' plastic elements for increment ',num2str(i), ...
         ' at load = ', num2str(i*dp)])
     es
   end;
   
   % Save variables for curve plotting
   pl(i+1,:)=[-a(6),i*dp];
 
end;
   % Plot force-displacement relation
   figure(1)
   plot(pl(:,1),pl(:,2),'-');
