function [ TimeDomainValues,LocationOfDataCarriers ] = OFDMModulation( FFTCoeff, FFTSize, DataCarriersCount )
%OFDMModulation Performs OFDMModulation of the Input Symbols
%   FFTCoeff -> Input Symbols after Digital Modulation.
%   FFTSize  -> Length of FFT
%   DataCarriesCount -> Number of useful data carriers.

%   Loads the Input Symbols into subcarriers with index
%   -DataCarriersCount/2 to DataCarriersCount/2 excluding DC
%   Then performs IFFT to obtain the OFDM time domain signal


%% Aid to load QPSK symbols in their corresponding sub carrier  %%
InputSequenceP = zeros(1,FFTSize/2);
InputSequenceN = zeros(1,FFTSize/2);
%% Loading QPSK symbols           %%
for(i=2:(DataCarriersCount/2)+1)
    InputSequenceP(i) = FFTCoeff(i-1);
end

for(i=1:(DataCarriersCount/2))
    InputSequenceN(i) = FFTCoeff(i+(DataCarriersCount/2));
end
%% Frequency Domain values of QFDM Signal %%
FFTCoeffLoaded = [ InputSequenceP flip(InputSequenceN)];
%% Time domain Signal of QFDM Signal      %%
TimeDomainValues = ifft(FFTCoeffLoaded,FFTSize)*sqrt(FFTSize);
LocationOfDataCarriers = find((FFTCoeffLoaded ~= 2) & (FFTCoeffLoaded ~= 0));

end

