function [f,delta] = bs_call(price,strike,int_rate,expiry,vol,size)
    lso = (log(price/strike)+(int_rate+(power(vol,2))/2)*expiry);
    d1 = lso./(vol*sqrt(expiry));
    d2 = d1 - vol*sqrt(expiry);
    f = size*(price.*normcdf(d1)-strike*exp(-int_rate*expiry)*normcdf(d2));
    delta = -size*normcdf(d1);
end