function [Xhead_up,Yhead_up,Zhead_up, Xhead_down,Yhead_down,Zhead_down] = Get_BF_IdealSperm_Head(scale_arc)


% Generate ideal head surface.
[Xhead_up,Yhead_up,Zhead_up, Xhead_down,Yhead_down,Zhead_down] = Generate_Virtual_Head_UpDownSides;  %xup1: M1*M2.


Xneck=max(max(Xhead_up)); 
Xhead_up=Xhead_up-Xneck;
Xhead_down=Xhead_down-Xneck;
  
 
Xhead_up = Xhead_up*scale_arc;
Yhead_up = Yhead_up*scale_arc; 
Zhead_up = Zhead_up*scale_arc;
Xhead_down = Xhead_down*scale_arc;  
Yhead_down = Yhead_down*scale_arc;  
Zhead_down = Zhead_down*scale_arc;

 