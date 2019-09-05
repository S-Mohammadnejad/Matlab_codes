function [ M, Models, num_models ] = Read_sequence( logfiles,num_logs,rxxmols,num_rxxmols )
%Read_sequence function: read a sequence of files (logfile and rxxmol)
%  put all logfiles and rxxmols under the same folder, e.g. 'MG_Ar_Pt'
%% Read logfiles
M=zeros(num_logs,1); %read logfiles to get the number of frames of each model
for i=1:num_logs
    [ M(i) ] = readlogfile_TAC( logfiles(i).name );
end

%% Read rxxmols
Models=cell(1,num_rxxmols); %a collection of all models' rxxmols
for i=1:num_rxxmols
    [ Models{i} ]=readrxxmol(rxxmols(i).name,M(i));
end

if num_rxxmols==num_logs %number of rxxmol files and logfiles should be the same
    num_models=num_rxxmols; 
end

end

function [ M ] = readlogfile_TAC( fname )
%readlogfile_TAC: Read log file
  
fid=fopen(fname);
if fid==-1
    error('File %s does not exist',fname);
end

parameters=textscan(fid, '%s',9,'Headerlines',1);
parameter=parameters{1};
value=fscanf(fid,'%f %f %f %f %f %f %f %f %*s',[length(parameter)-1,inf]);
value=value';

iteration=value(:,strcmpi('Iteration',parameter)); %iteration number
M=length(iteration); %number of frames

fclose(fid);

end

function [ Model ] = readrxxmol( fname,M )
%readrxxmol function: read rxxmol file.
%   fname: .rxxmol file name, must be in 'character string'.
fid=fopen(fname,'r');
if fid==-1
    error('File %s does not exist',fname);
end

%Model(i): a specific frame
for i=1:M
    MolNo=textscan(fid,'%f',1); %fetch molecule/atom number in one frame
    Model(i).MolNo=cell2mat(MolNo(1)); %numerical matrix
    a=textscan(fid,'%s %f %f %f %f %f %f %f',Model(i).MolNo,'Headerlines',2);
    Model(i).Element=a{1}; %character string cell
    Model(i).Coordinate=cell2mat(a(2:4)); %numerical matrix
    Model(i).Velocity=cell2mat(a(5:7)); %numerical matrix
    Model(i).Index=cell2mat(a(8)); %atom indices
end
fclose(fid);
end