fclose all;
close all;
clear;
clc;

norder =			0;
kind =				1;
nzeros =			10001;
filewrit = 			0;

roots =				besselzero(norder,nzeros,kind);

if (filewrit == 1)
	fid1 =			fopen(['besselroots_',num2str(kind),'_',num2str(norder),'_',num2str(nzeros),'.dat'],'w');
	fprintf(fid1,	[num2str(roots,'%.16d\t\t'),'\n']);
	fclose(fid1);
end

if (kind == 1)
	r =		linspace(0,5000,10001);
	x =		besselj(norder,r);
	figure();
	plot(r,x);
	hold on;
	plot(r,0*x);
	grid on;
end
