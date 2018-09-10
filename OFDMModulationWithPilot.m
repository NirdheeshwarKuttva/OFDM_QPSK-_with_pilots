function [ TimeDomainValues,LocationOfDataCarriers,PilotLocations ] = OFDMModulationWithPilot( FFTCoeff, FFTSize, DataCarriersCount )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
j = 2;
PilotSymbol = 5;
InputSequenceP = zeros(1,FFTSize/2);
InputSequenceN = zeros(1,FFTSize/2);
for(idx=1:FFTSize/2)
    if(rem(idx-1,9) == 0)
       InputSequenceP(idx) = PilotSymbol;
    else
        InputSequenceP(idx) = FFTCoeff(j-1);
        j=j+1;
        if(j>(DataCarriersCount/2)+1)
               break;
        end;
    end
end

j=1;        
for(idx=1:FFTSize/2)
    if(rem(idx-1,9) == 0)
       InputSequenceN(idx) = PilotSymbol;
    else
    InputSequenceN(idx) = FFTCoeff(j+(DataCarriersCount/2));
    j=j+1;
    if(j>(DataCarriersCount/2))
          break;
    end
    end
end
%% Frequency Domain values of QFDM Signal %%
FFTCoeffLoaded = [ InputSequenceP flip(InputSequenceN)];
FFTCoeffLoaded(FFTSize/2) = PilotSymbol;
LocationOfDataCarriers = find((FFTCoeffLoaded ~= 5) & (FFTCoeffLoaded ~= 0));
PilotLocations = find(FFTCoeffLoaded == 5);
VirtualCarriers = find(FFTCoeffLoaded == 0);
%% Time domain Signal of QFDM Signal      %%
TimeDomainValues = ifft(FFTCoeffLoaded)*sqrt(FFTSize);


end

