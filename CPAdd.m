function [ OFDMSymbol ] = CPAdd( TimeDomainValues,CPLen,BinLen )
%CPAdd Adds Cyclic Prefix
%    TimeDomainValues -> Signal to which CP needs to be added
%    CPLen-> Length of CP
%    FFTSize -> Size of FFT (techincally Number of Bins)
%    OFDMSymbol -> Signal with CP
%    Adds the last CPLen samples of Input signal as prefix to Input signal

OFDMSymbol = [TimeDomainValues((BinLen-CPLen+1):BinLen) TimeDomainValues];


end

