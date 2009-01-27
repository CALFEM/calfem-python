# example exs6 
#----------------------------------------------------------------
# PURPOSE 
#    Analysis of a plane frame.
#----------------------------------------------------------------

# REFERENCES
#     G"oran Sandberg 94-03-08 
#     Karl-Gunnar Olsson 95-09-28
#     Anders Olsson 99-03-01
#     Ola Dahlblom 2004-09-14
#----------------------------------------------------------------

from numpy import *

#----- Topology -------------------------------------------------

 Edof=[1  4  5  6 1  2  3;
       2  7  8  9 10 11 12;
       3  4  5  6  7  8  9];      

#----- Stiffness matrix K and load vector f ---------------------

 K=zeros(12);	f=zeros(12,1);	  f(4)=2e+3;

#----- Element stiffness and element load matrices  -------------

 E=200e9;
 A1=2e-3;     A2=6e-3;
 I1=1.6e-5;	  I2=5.4e-5;


 ep1=[E A1 I1];	 ep3=[E A2 I2];
 ex1=[0 0];      ex2=[6 6];	     ex3=[0 6];
 ey1=[4 0];      ey2=[4 0];	     ey3=[4 4];
 eq1=[0 0];
 eq2=[0 0];
 eq3=[0 -10e+3];

 Ke1=beam2e(ex1,ey1,ep1);
 Ke2=beam2e(ex2,ey2,ep1);
 [Ke3,fe3]=beam2e(ex3,ey3,ep3,eq3);

#----- Assemble Ke into K ---------------------------------------

 K=assem(Edof(1,:),K,Ke1);
 K=assem(Edof(2,:),K,Ke2);
 [K,f]=assem(Edof(3,:),K,Ke3,f,fe3);

#----- Solve the system of equations and compute reactions ------

 bc=[1 0; 2 0; 3 0; 10 0; 11 0];	
 [a,r]=solveq(K,f,bc)

#----- Section forces -------------------------------------------

 Ed=extract(Edof,a);

 es1=beam2s(ex1,ey1,ep1,Ed(1,:),eq1,21) 
 es2=beam2s(ex2,ey2,ep1,Ed(2,:),eq2,21) 
 es3=beam2s(ex3,ey3,ep3,Ed(3,:),eq3,21) 

 #----- Draw deformed frame ---------------------------------------
 
 figure(1)
 plotpar=[2 1 0];
 eldraw2(ex1,ey1,plotpar);
 eldraw2(ex2,ey2,plotpar);
 eldraw2(ex3,ey3,plotpar);
 sfac=scalfact2(ex3,ey3,Ed(3,:),0.1);
 plotpar=[1 2 1];
 eldisp2(ex1,ey1,Ed(1,:),plotpar,sfac);
 eldisp2(ex2,ey2,Ed(2,:),plotpar,sfac);
 eldisp2(ex3,ey3,Ed(3,:),plotpar,sfac);
 axis([-1.5 7.5 -0.5 5.5]); 
 pltscalb2(sfac,[1e-2 0.5 0]);
 axis([-1.5 7.5 -0.5 5.5]);
 title('displacements')
 
#----- Draw normal force diagram --------------------------------
 
 figure(2)
 plotpar=[2 1];
 sfac=scalfact2(ex1,ey1,es1(:,1),0.2);
 eldia2(ex1,ey1,es1(:,1),plotpar,sfac);
 eldia2(ex2,ey2,es2(:,1),plotpar,sfac);
 eldia2(ex3,ey3,es3(:,1),plotpar,sfac);
 axis([-1.5 7.5 -0.5 5.5]);
 pltscalb2(sfac,[3e4 1.5 0]);
 title('normal force')

#----- Draw shear force diagram ---------------------------------
 
 figure(3)
 plotpar=[2 1];
 sfac=scalfact2(ex3,ey3,es3(:,2),0.2);
 eldia2(ex1,ey1,es1(:,2),plotpar,sfac);
 eldia2(ex2,ey2,es2(:,2),plotpar,sfac);
 eldia2(ex3,ey3,es3(:,2),plotpar,sfac);
 axis([-1.5 7.5 -0.5 5.5]);
 pltscalb2(sfac,[3e4 0.5 0]);
 title('shear force') 

#----- Draw moment diagram --------------------------------------
 
 figure(4)
 plotpar=[2 1];
 sfac=scalfact2(ex3,ey3,es3(:,3),0.2);
 eldia2(ex1,ey1,es1(:,3),plotpar,sfac);
 eldia2(ex2,ey2,es2(:,3),plotpar,sfac);
 eldia2(ex3,ey3,es3(:,3),plotpar,sfac);
 axis([-1.5 7.5 -0.5 5.5]);
 pltscalb2(sfac,[3e4 0.5 0]);
 title('moment') 

#------------------------ end -----------------------------------
 echo off