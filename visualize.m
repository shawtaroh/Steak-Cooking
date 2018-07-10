function [x y] = visualize(h) %creates x and y coordinate matrices reflecting non-uniform h tensor
[Ny Nx d]=size(h);
x=zeros(Ny,Nx);
y=zeros(Ny,Nx);
for i=1:Nx
	if mod(Ny,2)==0 %if Ny is even, the origin should lie between the two center rows
    	y(Ny/2+1,i)=0.5*h(Ny/2+1,i,2);
        for j=Ny/2+2:Ny
			y(j,i)=y(j-1,i)+0.5*(h(j-1,i,2)+h(j,i,2));
		end
        y(Ny/2,i)=-0.5*h(Ny/2,i,2);
		for j=Ny/2-1:-1:1
			y(j,i)=y(j+1,i)-0.5*(h(j+1,i,2)+h(j,i,2));
		end
    else %if Ny is odd, the origin should lie on the center row
		for j=(Ny+1)/2+1:Ny
			y(j,i)=y(j-1,i)+0.5*(h(j-1,i,2)+h(j,i,2));
		end
		for j=(Ny+1)/2-1:-1:1
			y(j,i)=y(j+1,i)-0.5*(h(j+1,i,2)+h(j,i,2));
		end
    end
end
for j=1:Ny
	if mod(Nx,2)==0 %if Nx is even, the the origin should lie between the two center columns
    	x(j,Nx/2+1)=0.5*h(j,Nx/2+1,1);
		for i=Nx/2+2:Nx
			x(j,i)=x(j,i-1)+0.5*(h(j,i-1,1)+h(j,i,1));
		end
        x(j,Nx/2)=-0.5*h(j,Nx/2,1);
		for i=Nx/2-1:-1:1
			x(j,i)=x(j,i+1)-0.5*(h(j,i+1,1)+h(j,i,1));
		end
    else %if Nx is odd, the origin should lie on the center column
    	for i=(Nx+1)/2+1:Nx
			x(j,i)=x(j,i-1)+0.5*(h(j,i-1,1)+h(j,i,1));
		end
		for i=(Nx+1)/2-1:-1:1
			x(j,i)=x(j,i+1)-0.5*(h(j,i+1,1)+h(j,i,1));
		end
    end
end

%x(1,1)=(x(1,2)+x(2,1))/2;
%x(end,1)=(x(end,2)+x(end-1,1))/2;
%x(1,end)=(x(1,end-1)+x(2,end))/2;
%x(end,end)=(x(end,end-1)+x(end-1,end))/2;

%y(1,1)=(y(1,2)+y(2,1))/2;
%y(end,1)=(y(end,2)+y(end-1,1))/2;
%y(1,end)=(y(1,end-1)+y(2,end))/2;
%y(end,end)=(y(end,end-1)+y(end-1,end))/2;

end