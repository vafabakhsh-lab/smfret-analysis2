function [donor, acceptor] = readTraces(Ntraces, len, Data)

%{ 

This function ...

    Parameters
    ----------


    Returns
    -----------

%}





time=(0:(len-1))*timeunit;



% convert into donor and acceptor traces
index=(1:Ntraces*len);
Data=zeros(Ntraces,len);
donor=zeros(Ntraces/2,len);
acceptor=zeros(Ntraces/2,len);
Data(index)=raw(index);

for j=1:(Ntraces/2)
    donor(j,:)=Data(j*2-1,:);
    acceptor(j,:)=Data(j*2,:);
end
    

