%% QFDM-QPSK  

clc;
clear all;
close all;

%% Constants and Variables used throught .
Fs = 10.24*10^6; % Sampling Frequency Hz                   
FFTSize = 1024;  % FFT Size
Ts = 1/Fs;       % Sampling Period Seconds
T = 1/FFTSize;   % Symbol Time
OFDMSymbolCount=8;% Number of OFDM Symbols per Subframe
BitsPerSymbol= 2; % Number of bits per Symbol. 2 for QPSK
DataCarrierCount=900; % Number of Data carriers per OFDM symbol
TotalInputBitCount=BitsPerSymbol*DataCarrierCount* OFDMSymbolCount;%Total Number of Bits required to be generated
TotalSymbolCount  = TotalInputBitCount/BitsPerSymbol; % Total Number of Symbols
CPLen = 256;     % Length of Cyclic Prefix
PilotSymbol = 5;  % Value of the pilot symbol. 
OFDMWithPilot = [-1]; %OFDM symbols that carry Pilot symbols. Use -1 if no Training symbols are use.
PDPdb = [0,-3,-5,-6,-10];
PDP   = 10.^(PDPdb/10);
PDPNormalised = PDP/sum(PDP);

%% Generating Input Bits.
% rand function is used to generate TotalInputBitCount number of Inputs
% with values ranging from 0 to 1. Then it is rounded off to the nearest
% integer giving us bits.
BER = zeros(1,7);
for MC = 1:1      
InputData     = rand(1,TotalInputBitCount); % Generates Random Inputs 
RInputData    = round(InputData);           % Rounding Random Inputs    

%%  QPSK Modulation
% NRZ encoding is performed on the input bits followed by extracting
% Inphase and Quadrature components (odd and even respectively) to form
% QPSK symbols.

RInputData    = (2*RInputData - 1);         % NRZ encoding.
%  Variables used for QPSK Modulation 
SI = []; % Array to store the inphase bits
SQ = []; % Array to store the quadrature bits
FFTCoeff = []; % Array to store the QPSK symbols

SI = RInputData(1:2:TotalInputBitCount); % Odd bits are loaded onto the Inphase array
SQ = RInputData(2:2:TotalInputBitCount); % Even bits are loaded onto the Quadrature array

FFTCoeff = SI + 1i.*SQ;      % QPSK Symbols generation
FFTCoeff = reshape(FFTCoeff,[DataCarrierCount,OFDMSymbolCount]);%
FFTCoeffLoadedReshaped = zeros(FFTSize,OFDMSymbolCount); % Matrix to store the load FFTCoeff

%% OFDM Modulation
% Here the QPSK symbols are loaded onto the correct frequency bins.Also,
% Pilot symbols in OFDM symbol number 0 and 4 are added. Refer
% OFDMModulationWithPilot and OFDMModulation files for more information.
% Then Cyclic prefix of pre defined length is added.
% TimeDomainValues     -> Modulated time domain OFDM signal.
% DataCarriersLocation -> Index(Location) of Data carriers.
% PilotLocation        -> Location of Pilot symbols.

Subframe = []; % Array to store the Subframe

for SymbolCount = 1:OFDMSymbolCount  
    
    if( sum(SymbolCount == OFDMWithPilot+1)) % OFDM Symbols 0 and 4 contains Pilot symbols.
     % The corresponding symbols are passed to the function which packs the
     % symbols into the correct bins and adds pilot symbols then  returns
     % the OFDM modulated time domain signal among other things.
     [TimeDomainValues,DataCarriersLocation,PilotLocation] = OFDMModulationWithPilot(FFTCoeff ...
                        (:,SymbolCount),FFTSize,DataCarrierCount,PilotSymbol);
  
    else   
      % This function is similar to the previous one. Only difference is
      % that it does not load pilot symbols onto the frequency bins. This
      % function also returns the time domain OFDM Signal
     [TimeDomainValues,DataCarriersLocationNP] = OFDMModulation(FFTCoeff ...
                        (:,SymbolCount),FFTSize,DataCarrierCount);
    end
    
    FFTCoeffLoaded = fft(TimeDomainValues)/sqrt(FFTSize);
    
    % Cyclic prefix of length defined by variable CPLen is added.
    % Refer CPAdd.m file for more information 
    OFDMSymbol = CPAdd(TimeDomainValues,CPLen,FFTSize);
    
    FFTCoeffLoadedReshaped(:,OFDMSymbolCount) = FFTCoeffLoaded;
    
    Subframe = [Subframe OFDMSymbol];% The Time domain OFDMSymbols are concatenated to form the Subframe array
end

OFDMSymLen = FFTSize+CPLen; % OFDM Symbol length
SubframeLen = (OFDMSymLen)*OFDMSymbolCount;% Subframe length
subplot(3,2,1)
plot((0:SubframeLen-1)*Ts,real(Subframe));
title(' OFDM Subframe ');
xlabel(' Time (Seconds)');
ylabel(' Amplitude ');

%% BER Calculations
% Here the Bit error rate vs SNR  of OFDM in a AWGN channel is calculated and plotted.
% RecievedSubframe  -> Subframe after probagating through the channel.
% OFDMSymbol        -> Contains the extracted Recieved OFDM symbols
% OSWCP             -> Contains the OFDM data points that will be fed to the FFT block
% UsefulQPSKSymbols -> Contains the Recieved Modulated Data symbols
% RecievedSI        -> Recieved Inphase symbols
% RecievedSQ        -> Recieved Quadrature symbols

SNRdb = 0:1:12; % SNR range.
SNRLinear  = 10.^(SNRdb/10);% Eb/No
RecievedFFTCoeff = zeros(1,DataCarrierCount*SymbolCount);% Array to store the recieved symbols
plotIdx = 2; % Subplot number 
Index   = 1;
for k=-2:2; % Synchronisation error. 
    if( OFDMWithPilot(1) == -1)
        H = exp(1i*2*pi*[1:FFTSize]*k/FFTSize);
        HInv = 1./H;
    end
    for snr = SNRdb;
        RecievedSubframe = awgn(Subframe,snr); % AWGN is added to the signal.
         for l = 0:1:OFDMSymbolCount-1
            % Here OFDM symbols are extracted from the subframe one by one
            % and then the Synchronisation error is mimicked 
            
            OFDMSymbol =  RecievedSubframe(l*(OFDMSymLen)+1: (l+1)*OFDMSymLen);% Extracting OFDM symbols from Subframe           
            if(k <= 0)  
                OSWCP      = OFDMSymbol(CPLen+1+k:OFDMSymLen+k); % Synchronisation Error. K=0 => No error
            else % k > 0
                % For certain cases, Some points from the next OFDM symbol
                % is appended to the end of the current FFT window to mimic
                % ISI. And in the case of last OFDM symbol, Zeros are
                % added.
                
                if(l ~= OFDMSymbolCount-1)
                     OSWCP      =[OFDMSymbol(CPLen+1+k:OFDMSymLen) RecievedSubframe((l+1)*(OFDMSymLen)+1: (l+1)*(OFDMSymLen)+k)];
                else
                     OSWCP      =[OFDMSymbol(CPLen+1+k:OFDMSymLen) zeros(1,k)];
                end
            end
            
           
            
            QPSKSymbols =( fft(OSWCP,FFTSize)/sqrt(FFTSize));% Finding the Frequency domain values
            
            % Channel estimation using Pilot Symbols
             if(sum(l == OFDMWithPilot) && OFDMWithPilot(1) ~= -1)
                   H = zeros(1,FFTSize);
                  H(PilotLocation) = (QPSKSymbols(PilotLocation)/PilotSymbol);
                    
                   H(find(H==0)) = interp1(PilotLocation,H(PilotLocation),find(H==0));
                   HInv = 1./H;
                    
%                     HInv = smooth(HInv,128);
%                     HInv = reshape(HInv,[1,FFTSize]);
            end
            
            QPSKSymbols = HInv.*QPSKSymbols ./ (abs(HInv));
            % Pilots and Virtual carriers are discarded.
            if(sum(l == OFDMWithPilot)&& OFDMWithPilot(1) ~= -1)
                UsefulQPSKSymbols = QPSKSymbols(DataCarriersLocation);
            else
                UsefulQPSKSymbols = QPSKSymbols(DataCarriersLocationNP);
            end
            % The latter half of the Modulated symbols are flipped to
            % compensate the FFTShift.
            UsefulQPSKSymbols= [UsefulQPSKSymbols(1:DataCarrierCount/2) flip(UsefulQPSKSymbols(1+DataCarrierCount/2:DataCarrierCount))];
            
            % The Recieved Modulated symbols are accumulated for ber
            % calculations.
            RecievedFFTCoeff((l*DataCarrierCount)+1:(l+1)*DataCarrierCount) = UsefulQPSKSymbols;
        end
        
        % Code for Equalisation here
        % Resphape RecievedFFTCoeff
        % Estimate the channel using first and fifth row.
       
        
        % QPSK Demoulation
        RecievedSI = sign(real(RecievedFFTCoeff)); % Extracting the Inphase components
        RecievedSQ = sign(imag(RecievedFFTCoeff)); % Extracting the Quadrature components
       
        ErrorCount = sum( (RecievedSI ~= SI) + (RecievedSQ ~= SQ)); % Total number of errors
        BER(Index) =  ErrorCount/TotalInputBitCount; % Bit error rate
        Index = Index+1;
    end
    Index = 1;
        subplot(3,2,plotIdx)
    BER = BER ;
    tber = 0.5*erfc(sqrt(SNRLinear)); % Theorectial Bit error rate of QPSK
     
    semilogy(SNRdb,BER,'-bo',SNRdb,tber,'-mh'); % Plotting the theorectical and simulated BER
    xlabel(' SNR(db)');
    ylabel(' BER ');
    t = ['With Sync Error k = ',num2str(k)];
    title(t);
    plotIdx = plotIdx+1;
end
   
end

    


        
  