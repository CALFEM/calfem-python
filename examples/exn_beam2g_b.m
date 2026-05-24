% example exn_beam2g_b
%----------------------------------------------------------------
% PURPOSE 
%    Buckling analysis of a plane frame.
%----------------------------------------------------------------

% LAST MODIFIED: O Dahlblom  2015-12-04
% Copyright (c)  Division of Structural Mechanics and
%                Division of Solid Mechanics.
%                Lund University
%--------------------------------------------------------------------------
 echo off
 
% ----- Topology ----------------------------------------------------------

 Edof=[1  4  5  6  1  2  3;
       2  7  8  9 10 11 12;
       3  4  5  6  7  8  9]; 

% ----- Element properties and global coordinates ------------------------- 
      
 E=200e9;  
 A1=2e-3;       A2=2e-3;       A3=6e-3;
 I1=1.6e-5;     I2=1.6e-5;     I3=5.4e-5;
 ep1=[E A1 I1];	ep2=[E A2 I2]; ep3=[E A3 I3];   
 eq1=[0];       eq2=[0];       eq3=[-50e3];
 ex1=[0 0];     ex2=[6 6];     ex3=[0 6];   
 ey1=[4 0];     ey2=[4 0];     ey3=[4 4];

% ----- Initial axial forces ----------------------------------------------

   QX1=0.0001;   QX2=0;        QX3=0;
   QX01=1;

% ----- Iteration for convergence -----------------------------------------

   eps=1e-6;		
   n=0;			
   while(abs((QX1-QX01)/QX01)>eps)
     n=n+1;

     K=zeros(12,12);
     f=zeros(12,1);	
     f(4)=10e3;

     [Ke1]=beam2ge(ex1,ey1,ep1,QX1);
     [Ke2]=beam2ge(ex2,ey2,ep2,QX2);
     [Ke3,fe3]=beam2ge(ex3,ey3,ep3,QX3,eq3);
     
     K=assem(Edof(1,:),K,Ke1);
     K=assem(Edof(2,:),K,Ke2);
     [K,f]=assem(Edof(3,:),K,Ke3,f,fe3);
     if n==1;
        K0=K;
     end;    
    
     bc=[1 0;2 0;3 0;10 0;11 0];	
     [a,r]=solveq(K,f,bc);

     Ed=extract_ed(Edof,a);

     QX01=QX1; 
     [es1,QX1,edi1]=beam2gs(ex1,ey1,ep1,Ed(1,:),QX1,eq1,21);
     [es2,QX2,edi2]=beam2gs(ex2,ey2,ep2,Ed(2,:),QX2,eq2,21);
     [es3,QX3,edi3]=beam2gs(ex3,ey3,ep3,Ed(3,:),QX3,eq3,21);
     
     if(n>20)
        disp('The solution does not converge')
        break
     end
   end

% ----- Buckling analysis -------------------------------------------------
 echo on  
 
 b=bc(:,1);
 [lambda,phi]=eigen(K,K0,b);
 nmods=size(lambda);
 one=ones(nmods);
 alpha=one./(one-lambda);
 alpha(1)
 phi(:,1)
  
 %----- Draw shape at instability -----------------------------------------
 
 Ed=extract_ed(Edof,-phi(:,1));
 [es1,QX1,edi1]=beam2gs(ex1,ey1,ep1,Ed(1,:),[0],[0],21)
 [es2,QX2,edi2]=beam2gs(ex2,ey2,ep2,Ed(2,:),[0],[0],21)
 [es3,QX3,edi3]=beam2gs(ex3,ey3,ep3,Ed(3,:),[0],[0],21)
 figure(1)
 plotpar=[3 1 0];
 eldraw2(ex1,ey1,plotpar);
 eldraw2(ex2,ey2,plotpar);
 eldraw2(ex3,ey3,plotpar);
 sfac=scalfact2(ex3,ey3,edi3,0.1);
 plotpar=[1 2 0];
 dispbeam2(ex1,ey1,edi1,plotpar,sfac);
 dispbeam2(ex2,ey2,edi2,plotpar,sfac);
 dispbeam2(ex3,ey3,edi3,plotpar,sfac);
 axis([-1.5 7.5 -0.5 5.5]); 
 title('Shape at instability')
 
 figure(2)
 sfac=scalfact2(ex3,ey3,Ed(3,:),0.1);
 plotpar=[3 1 0];
 eldraw2(ex1,ey1,plotpar);
 eldraw2(ex2,ey2,plotpar);
 eldraw2(ex3,ey3,plotpar);
 Ed=extract_ed(Edof,-phi(:,1));
 plotpar=[1 2 0];
 eldisp2(ex1,ey1,Ed(1,:),plotpar,sfac);
 eldisp2(ex2,ey2,Ed(2,:),plotpar,sfac);
 eldisp2(ex3,ey3,Ed(3,:),plotpar,sfac);
 axis([-1.5 7.5 -0.5 5.5]); 
 title('Shape at instability 1')
 
 echo off

