function [z] = impute(y1, y2)
    n =  length(y1);
    z =  zeros(n,1);
    for i = 1:n
        if (y1(i)>0) && (y2(i)>0)
            z(i) = (y1(i)+y2(i))/2;
        elseif (y1(i)>0) 
            z(i) = y1(i);
        elseif (y2(i)>0)
            z(i) = y2(i);
        else
            z(i) = -99;
        end
    end
    z = z(z>0);
end