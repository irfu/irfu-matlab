function res = find_zero(y,i1,i2)

a = ( y(i1) - y(i2) ) / ( i1 - i2 );
b = ( y(i1) + y(i2) - a*( i1 + i2 ) ) / 2;

res = -b/a;
end