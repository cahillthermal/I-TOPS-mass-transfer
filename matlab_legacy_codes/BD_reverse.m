function [ f_inver,data_inver ] = BD_reverse( f_collect,data_collect )
%This reverse the data in frequency (used when the collecting data from high frequency)
f_inver=f_collect.*0;
data_inver=data_collect.*0;

Nmax=length(f_collect);
for i=1:1:Nmax
f_inver(i) = f_collect(Nmax+1-i);
data_inver(i) = data_collect(Nmax+1-i);
end

