fclose all;
close all;
clear;
clc;

type =	'Y';
N =		10001;

kernel = zeros(N,N);

if (strcmp(type,'Y'))
	for ix = 1:N
		for ir = 1:N
			kernel(ix,ir) =		(2.0/(
		end
	end
elseif (strcmp(type,'K'))
	
else
	disp('Wrong Kernel type supplied, not calculating and writing any kernel');	
end
