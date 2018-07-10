meanPhi = 0;
area = 0;
for i = 1:P.Nx
    for j = 1:P.Ny
        meanPhi = meanPhi + S.phi(j,i)*h(j,i,1)*h(j,i,2);
        area = area + h(j,i,1)*h(j,i,2);
    end 
end

meanPhi = meanPhi/area;
meanN = (1-meanPhi)/(1+.3*meanPhi)