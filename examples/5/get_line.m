function [line_str] = get_line( fname )

%fname= 'cmu.truth';
fptr = fopen(fname, 'r');
lines = 1;
attri = {};

while ~feof(fptr)
attri{lines} = fgetl(fptr);
lines = lines + 1;
end

fclose(fptr);

line_str = attri;

end