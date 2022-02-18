%{ 

    About

    Parameters
    ----------

    Returns
    ----------


    Changelog
    ----------

%}


%{ 

    Pseudocoding
    
    1. clean workspace
    2. declare variables
    3. get data directory
    4. get data parameters
    5. read data
    6. for each trace in data 
        a. calculate FRET
        b. graph trace
        c. ask to save or skip
    7. save trace
        a. select part of trace to save
        b. define background
        c. write data to file
    

%}


% close any open windows and clear any variables
close all;
fclose('all');
clearvars;

% input file directory and file name to be analyzed.
pth=input('Directory: ');
cd(pth);

% Specify how many points to average.
sm=3;
% Specify the Gamma factor
gamma=1.0;

gap=0:0.01:1;
sg=size(gap);sg=sg(2);
N_His=zeros(1,sg);
N_His_raw=zeros(1,sg);

fname=num2str(input('index # of filename [default=1]  '));
	if isempty(fname)
   	fname='1';
	end
fname=num2str(fname);
filename=['film' fname '.traces']
fid=fopen(['film' fname '.traces'],'r');

A=dir;
[nf,dumn]=size(A);
numberoffiles=length(A)-2;
disp(['Number of files: ',num2str(numberoffiles)]);

timeunit=input('Time unit: [default=0.06 sec] ');
    if isempty(timeunit)
        timeunit=0.06;
    end
    
    
leakage=input('Donor leakage correction: [default=0.09] ');
    if isempty(leakage)
      leakage=0.09;
    end
    
    
bg_definition_input=input('Manually define background? y or n [default=y] ');
    if isempty(bg_definition_input)
        bg_definition_input='y';
    end
    
bg_definition=1;

if bg_definition_input=='n'
    bg_definition_input=0;
end

saved_frames=zeros(1,100000);
number_of_frames=0;
total_N_of_frames=0;
was_deleted=0;


% repeat over all the files in the directory


end_of_directory=0;

    len=fread(fid,1,'int32'); 
    disp(['The len of the time traces is ',num2str(len)]);
    Ntraces=fread(fid,1,'int16');
    disp(['The number of traces is ',num2str(Ntraces/2)]);
   
    time=(0:(len-1))*timeunit;
    
    raw=fread(fid,Ntraces*len,'int16');
    fclose(fid);
    
%[donor, acceptor] = readTraces(Ntraces, len, timeunit, Data);
