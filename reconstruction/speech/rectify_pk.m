function [out] = rectify_pk(in)
    in(in<0) = 0;
    out = in;
end