function [Fre, Amp, Pha, PS] = Get_Vector_Calibrated_FFT(Function, nmode, Fs)

% This subfunction is to calculate FFT results (cumulative results in 'ns'...
% direction, by default) of the first 'nmode' dynamic FFT terms.
%
% Input: 'Function', ns*nt.
%        'nmode', 1*1, how many FFT dynamic modes are required, equal to 1 or 2.
%        'Fs', 1*1, sampling frequency, supposed to equal to 90 as our video framerate is 90fps.
%        
% Output: 'Fre', 1*nmode, frequency for the first 'nmode' dynamic FFT modes.
%         'Amp', ns*(nmode+1), amplitude for the static and the first 'nmode' dynamic FFT modes.
%         'Pha', ns*nmode, phase for the first 'nmode' dynamic FFT modes.
%         'PS', 1*(nt/2+1), cummulative power spectrum in the frequency
%         domain ( 'FD', 1* (round(nt/2)+1) ), not sorted.

% X_r(s,t) = X_fft0 +X_fft1 +X_fft2 +... 
% X_fft = Amp*cos(2*pi*Fre*t+Pha)


ns = size(Function,1);
nt = size(Function,2);
Amp = nan(ns,nmode+1);
Fre = nan(1,nmode);
Pha = nan(ns,nmode);
PS = nan(1,round(nt/2)+1);



%L=2^nextpow2(nt);  % Length of signal
L = nt;  % Length of signal
FD=Fs*(0:round(L/2))/L;  %Frequency domain of the FFT result, 1*(L/2+1) array.
F_DFT = fft(Function,L,2); % direct fourier transform of 'Function'. size(F_DFT)=[ns L]=[ns nt] 
P2_F = abs(F_DFT/L); % FFT power spectrum (calculate the two-sided spectrum P2). [ns*L] matrix
P1_F = P2_F(:,1:round(L/2)+1);  % calculate the one-sided spectrum P1 based on P2. [ns*(L/2+1)] matrix
P1_F(:,2:end-1) = 2*P1_F(:,2:end-1);  % calculate the one-sided spectrum P1 based on the even-numbered signal length L.



%% Preliminary sorting of the FFT modes, according to raw power spectrum 
% result, which however contains much noise and will thereby lead to inaccurate sorting.

%P1_F_cumsum_temp = cumsum(P1_F); % cummulative power along 'ns' direction, [ns*(L/2+1)] matrix.
%PS = P1_F_cumsum_temp(end,:); % cummulative power along 'ns' direction, [1*(L/2+1)].
PS = P1_F(end,:);   %//////////////////////////////////////////////////////////////////////////////////////
    %
    % I tried to smooth the power spectrum data to get the approariate 1st
    % and 2nd harmonic peaks (frequencies)  by using MATLAB function
    % 'smoothdata', but the smoothed power spectrum and resultant peaks barely 
    % have any difference. Therefore, I give up the smoothing.
    %Pcum = smoothdata(Pcum);
    %
[~, Ind_PS_sort] = sort(PS(2:end-1),'descend'); 
Ind_PS_sort = [1 Ind_PS_sort+1 round(nt/2)+1]; % index for the descending-sorted cummulative power, [1*(L/2+1)]. Noted that the first element is for the 0th FFT mode.
FD_sort = FD(Ind_PS_sort); % 1*(L/2+1)
 


%% Calibration step1: decide the 1st frequency (temporary) according to the 
% power spectrum, after removing the noise.

% The primary frequencies calculated here might have problems
% -- there may be power-spectrum peaks appearing very close to 0HZ, and this peak
% would be recognized by my code as the 1st FFT frequency, which is, however, 
% inappropriate. 

Ind_Fre1st_temp = find(FD_sort>1,1);
Ind_Fre1st = find(FD == FD_sort(Ind_Fre1st_temp));
Fre1st = FD(Ind_Fre1st);
%PS1st = PS(Ind_Fre1st);



if nmode==1  
    
%% Save data.

Fre = [Fre1st];


Pha_temp = angle(F_DFT);
Pha = Pha_temp(:,[Ind_Fre1st]);


Amp_temp = P1_F(:,[Ind_Fre1st]);
Amp_0mode = P1_F(:,1).*cos(Pha_temp(:,1));
Amp = [Amp_0mode Amp_temp];





elseif nmode==2

%% Calibration step2: decide the 2nd frequency (temporary) according to the
% power spectrum, which is supposed to be half or twice of the 1st frequency value.


% Locate the ideal 2nd peak of the power spectrum.

% Check if 'Fre1st_temp/2 ' is the 2nd peak.
    k1 = find(  abs( FD_sort-Fre1st/2 )<=1 );
    if isempty(k1)
        % 'Fre1st_temp/2' is not the 2nd peak we're looking for.
        PS2nd_V1 = 0;    
    else
        % 'Fre1st_temp/2' could be the 2nd peak we're looking for.  
        for i_k = 1:length(k1)
            Ind1_region(i_k) = find(FD==FD_sort(k1(i_k)));
        end
        [PS2nd_V1, ~]= max( PS(Ind1_region) );
        Ind_Fre2nd_V1 = find( PS==PS2nd_V1,1 );
    end

    
% Check if 'Fre1st_temp*2 ' is the 2nd peak.
    k2 = find(  abs( FD_sort-Fre1st*2 )<=1 );
    if isempty(k2)
        % 'Fre1st_temp*2' is not the 2nd peak we're looking for.
        PS2nd_V2 = 0;       
    else       
        % 'Fre1st_temp*2' could be the 2nd peak we're looking for.        
         for i_k = 1:length(k2)
            Ind2_region(i_k) = find(FD==FD_sort(k2(i_k)));
         end
         [PS2nd_V2, ~]= max( PS(Ind2_region) );
         Ind_Fre2nd_V2 = find( PS==PS2nd_V2,1 );
    end
    
    
% Determine which is the 2nd peak, 'Fre1st_temp/2' or 'Fre1st_temp*2'?
    if PS2nd_V1>=PS2nd_V2
        %PS2nd = PS2nd_V1;
        Fre2nd = FD(Ind_Fre2nd_V1);
    else
        %PS2nd = PS2nd_V2;
        Fre2nd = FD(Ind_Fre2nd_V2);
    end
    
    
% If without calibration for the 2nd frequency..   
   if 0
    Ind_Fre2nd_temp = find(FD_sort>1,2);  Ind_Fre2nd_temp = Ind_Fre2nd_temp(2);
    Ind_Fre2nd = find(FD == FD_sort(Ind_Fre2nd_temp));
    Fre2nd = FD(Ind_Fre2nd);
    PS2nd = PS(Ind_Fre2nd);
   end
    

   
   
%% Calibration step3: decide the (final) 1st and 2nd frequencies according to their frequency magnitudes.

% The final 1st frequency should be smaller than, and about half of, the
% 2nd frequency, while the results above may not satisfy this condition,
% thus mixing the 1st and 2nd frequencies.

if Fre1st<Fre2nd
    Fre1 = Fre1st;
    Fre2 = Fre2nd;
    %PS1 = PS1st;
    %PS2 = PS2nd;
else
    Fre1 = Fre2nd;
    Fre2 = Fre1st;
    %PS1 = PS2nd;
    %PS2 = PS1st;
end


Ind_Fre1st = find(FD==Fre1);
Ind_Fre2nd = find(FD==Fre2);



%% Save data.

Fre = [Fre1 Fre2];


Pha_temp = angle(F_DFT);
Pha = Pha_temp(:,[Ind_Fre1st Ind_Fre2nd]);


Amp_temp = P1_F(:,[Ind_Fre1st Ind_Fre2nd]);
Amp_0mode = P1_F(:,1).*cos(Pha_temp(:,1));
Amp = [Amp_0mode Amp_temp];


end