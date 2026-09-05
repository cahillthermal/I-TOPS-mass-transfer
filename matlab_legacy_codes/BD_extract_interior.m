function [f_interior,data_interior] = BD_extract_interior(f,data,f_min,f_max)
%EXTRACTS all points in range tmin<= t_data <= tmax, and corresponding
%points of f (of same length)

if f(1)>f(end)
[f,data] = BD_reverse(f,data); % reverse the data(when scanning from high frequency)
end

    Ndata=length(f);
    Nmin=Ndata;
    for i=Ndata:-1:1
        if f(i)>=f_min
            Nmin=i;
        end
    end
    Nmax=Nmin;
    for i=Nmin:Ndata
        if f(i)<=f_max
            Nmax=i;
        end
    end
    f_interior=f(1,Nmin:Nmax);
    data_interior=data(1,Nmin:Nmax);
end

