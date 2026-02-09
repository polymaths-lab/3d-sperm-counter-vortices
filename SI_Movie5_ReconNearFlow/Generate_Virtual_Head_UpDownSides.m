function [xup1,xup2,xup3, xdown1,xdown2,xdown3] = Generate_Virtual_Head_UpDownSides


% Sperm head:up+down£¬different colors,surf
    model.a1=2.0/45;
    model.a2=1.6/45;
    model.a3=1.0/45;
    nth=20;%10;
    nphi=40;%20;
    a=1;
    %M=nth*nphi;
    
    %head: up side
    [xup1,xup2,xup3]=GenerateSphereSurfaceForVisualisation_up(nth,nphi,a);
    %[M1 M2]=size(xup1);
    xup1=xup1*model.a1;
    xup2=xup2*model.a2;
    xup3=xup3*model.a3;
    
    %head: down side
    [xdown1,xdown2,xdown3]=GenerateSphereSurfaceForVisualisation_down(nth,nphi,a);
    xdown1=xdown1*model.a1;
    xdown2=xdown2*model.a2;
    xdown3=xdown3*model.a3;