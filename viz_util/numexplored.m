% to compute the rate at which the flyers explore the network

numnum = zeros(1,100);
for qqq=1:100
    [uniqpath,numuniqpath] = count_unique(reshape(flightpath(1:5,1:qqq),1,5*qqq));
    numnum(qqq) = size(numuniqpath,1);
end