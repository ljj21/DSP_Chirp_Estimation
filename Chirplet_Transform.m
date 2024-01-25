function [Spec,Freq,t] = Chirplet_Transform(Sig,fLevel,WinLen,SampFreq,alpha)
% This function is used to calculate the Chirplet transform of the signal.
% --------------------INPUT------------------%
% Sig    £ºOne-dimensional signal sequence to be analyzed
% fLevel £ºfrequency axis points associated with the Spec(in Bins)
% WinLen £ºGauss window width(in Bins)
% SampFreq£ºSignal sampling frequency(in Hz)
% alpha  £ºChirplet's line frequency modulation rate(in Hz/s)
% --------------------OUTPUT------------------%
% Spec   £º2D spectrum results (horizontal time axis, vertical frequency axis)
% Freq   £ºvertical frequency axis(in Hz)
% t      £ºhorizontal time axis(in Sec)
%--------------------------------------------------------%
% written by Guowei Tu, 28/07/2019 in SME,SJTU (Contact me via GuoweiTu@sjtu.edu.cn)
%--------------------------------------------------------%
% the Original Chirplet Transform Algorithm is introduced in [1];
% while the codes here is motivated by [2], which provides a new insight to the Chirplet Transform.
%--------------------------------------------------------%
% [1].Steve Mann, Simon Haykin, The Chirplet Transform: Physical Considerations, 
% IEEE TRANSACTIONS ON SIGNAL PROCESSING,VOL. 43, NO. 11, NOVEMBER 1995
% [2].Peng Z.K , Meng G., Lang Z.Q.,Chu F.L, Zhang W.M., Yang Y., Polynomial Chirplet Transform with Application to Instantaneous Frequency Estimation,
% IEEE Transactions on Measurement and Instrumentation 60(2011) 3222-3229

if (nargin < 1)
    error('At least one parameter required!');
end
if (nargin < 4)       % SampFreq default 1000Hz
    SampFreq = 1000;
end
if (nargin < 3)       % WinLen default 64 points
    WinLen = 64;
end
if (nargin < 2)       % fLevel defaults to 512 points
    fLevel = 512;
end

%% data preparation
SigLen = length(Sig); % Signal length
t = (0:1:SigLen-1)*(1/SampFreq); % time axis associated with the Signal
Sig = hilbert(Sig); % Calculate the analytical signal
%% Frequency axis and its half-axis points
fLevel = ceil(fLevel/2) * 2+1;    % Round the frequency axis length value fLevel in +¡Þ direction to the nearest odd number
Lf = (fLevel - 1)/2;    % Half frequency axis data points (fLevel has been rounded to an odd number)
%% Generate Gauss window functions
WinLen = ceil(WinLen/2) * 2+1;    % Round the window length value WinLen to the +¡Þ direction to the nearest odd number
WinFun = exp(-6* linspace(-1,1,WinLen).^2);    % Gauss window function, fixed time width [-1, 1], divided into WinLen modulation data points
Lw = (WinLen - 1)/2;    % Half window data points
%% CT spectrum result array (initialized to 0)
Spec = zeros(fLevel,SigLen); 
%% Sliding window signal,data segmentation preparation
for iLoop = 1:SigLen
    % Determine the upper and lower limits of the left and right signal index subscripts (note that to prevent the edge width from exceeding the time domain, the retrieval error!)
    iLeft = min([iLoop-1, Lw, Lf]);
    iRight = min([SigLen-iLoop, Lw, Lf]); 
    iIndex = -iLeft:iRight;
    
    iIndex1 = iIndex + iLoop;   % subscript index of the original signal
    iIndex2 = iIndex + Lw + 1;  % subscript index of window function vector
    Index = iIndex + Lf +1;     % Subscript index of the frequency axis (row number) in the Spec two-dimensional array
    
    R_Operator = exp(-1i * 2*pi*alpha * t(iIndex1).^2 / 2); % Frequency rotation operator (Shear Frequency)
    S_Operator = exp(1i * 2*pi*alpha * t(iLoop) * t(iIndex1)); % Frequency shift operator (Shift Frequency)
    
    Sig_Section = Sig(iIndex1).*R_Operator.*S_Operator; % Signal segment after frequency rotation and frequency shift
    
    Spec(Index, iLoop) = Sig_Section.* conj(WinFun(iIndex2));  % fill the two-dimensional array 
end
%% STFT
Spec = fft(Spec); % STFT
Spec = Spec*2/fLevel; % restores the true magnitude of FT
Spec = Spec(1:(end-1)/2,:); % till the Nyquist frequency

fk = 0:1:fLevel-1;
fk = fk(1:(end-1)/2);
Freq = linspace(0,0.5*SampFreq,length(fk)); % Output frequency axis for plotting
end

