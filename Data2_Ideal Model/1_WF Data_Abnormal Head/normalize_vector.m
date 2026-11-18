% Normalizing vector  

function [Tx,Ty,Tz] = normalize_vector(tx,ty,varargin)


if length(varargin)>=1
    tz  = varargin{1};
    
    t_modulus = (tx.^2 +ty.^2 +tz.^2).^(0.5);
    Tx        = tx./t_modulus;
    Ty        = ty./t_modulus;
    Tz        = tz./t_modulus;
    
else
    t_modulus = (tx.^2 +ty.^2).^(0.5);
    Tx        = tx./t_modulus;
    Ty        = ty./t_modulus;
end




% Old function

if 0 %<-------------------------|
t_modulus = (tx.^2 +ty.^2 +tz.^2).^(0.5);

Tx=tx./t_modulus;
Ty=ty./t_modulus;
Tz=tz./t_modulus;

end  %<-------------------------|