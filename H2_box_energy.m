% Computing different types of energy components for H2 gas
clear
close all
clc
%% Input parameters
kB=1.38064852*10^(-23); %[SI unit]
m_H2=2*1.0079/1000/(6.02*10^23); %molecular mass [kg]
mg=m_H2; %<<decide which TAC to be calculated: argon or hydrogen
mu=m_H2/4; %reduced mass for H2 molecule
%<<change temperatures
Tg=295; %gas molecule shot/initial temperature [K]
%v_2_w=4*kB*Tw/mg;
%v_mp_Tw=sqrt(2*kB*Tw/mg);
%% Read files
logfiles=dir('*.logfile');
num_logs=length(logfiles);
rxxmols=dir('*.rxxmol');
num_rxxmols=length(rxxmols);
[ num_frame, Models, num_try ] = Read_sequence( logfiles,num_logs,rxxmols,num_rxxmols ); 
fname=logfiles(1).name;
N_atom_gas=Models{1,1}(1).MolNo;
T=compute_Tm(fname);
T_mean=mean(T(:,2));
E_g_tr=1.5*kB*T_mean;
E_g_ro=kB*T_mean;
E_g_tot=2.5*kB*T_mean;
%% Extracting velocity and coordinate of gas atom
H2_data=zeros(num_try*N_atom_gas*num_frame,6); %multiply by 2: Hydrogen has 2 atoms
counter=1;
cc=0;
for i=1:num_try
    for j=1:num_frame(1)
        
            H2_data(counter:(cc+N_atom_gas),1:3)=Models{1,i}(j).Coordinate(1:N_atom_gas,:);
            H2_data(counter:(cc+N_atom_gas),4:6)=Models{1,i}(j).Velocity(1:N_atom_gas,:);
            cc=cc+N_atom_gas;
            counter=cc+1;
     
    end
end
%% Computeing C.O.M properties for H2
H2_com=zeros(size(H2_data));
for j=1:2:length(H2_data)
    H2_com(j,1)=(H2_data(j,1)+H2_data(j+1,1))/2;
    H2_com(j,2)=(H2_data(j,2)+H2_data(j+1,2))/2;
    H2_com(j,3)=(H2_data(j,3)+H2_data(j+1,3))/2;
    H2_com(j,4)=(H2_data(j,4)+H2_data(j+1,4))/2;
    H2_com(j,5)=(H2_data(j,5)+H2_data(j+1,5))/2;
    H2_com(j,6)=(H2_data(j,6)+H2_data(j+1,6))/2;
end
n=find(H2_com(:,1));
H2_com=H2_com(n,:);
%% computing velocity components for angular velocity
 rot_v=zeros(size(H2_data,1),12);
 for j=1:2:length(H2_data)
     [v_t11,v_t12,v_t21,v_t22]=Angular_Velocity(H2_data(j,1:3),H2_data(j,4:6),H2_data(j+1,1:3),H2_data(j+1,4:6));
     rot_v(j,:)=[v_t11,v_t12,v_t21,v_t22];
          
 end
 n=find(rot_v(:,1));
rot_v=rot_v(n,:);
 %% Rotational energy
 energy_rot=zeros(size(rot_v,1),1); %for each molecule
 for i=1:size(rot_v,1)
     energy_rot(i)=0.5*mu*(norm(rot_v(i,1:3)-rot_v(i,7:9)).^2+norm(rot_v(i,4:6)-rot_v(i,10:12)).^2);
 end
 %%%%%%%%%%%%%%%%%%%%%%%%%
 energy_rot_frame=zeros(num_frame,1);
 rr=length(energy_rot)/num_frame;
 cc=0;
 counter=1;
 for i=1:num_frame
     energy_rot_frame(i)=mean(energy_rot(counter:cc+rr));
     cc=cc+rr;
     counter=cc+1;
 end
 E_rot_m=mean(energy_rot_frame);
  %% Translational energy
 energy_tr=zeros(size(H2_com,1),1);
  for i=1:size(H2_com,1)
     energy_tr(i)=0.5*m_H2*norm(H2_com(i,4:6))^2;  
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  energy_tr_frame=zeros(num_frame,1);
 rr=length(energy_tr)/num_frame;
 cc=0;
 counter=1;
 for i=1:num_frame
     energy_tr_frame(i)=mean(energy_tr(counter:cc+rr));
     cc=cc+rr;
     counter=cc+1;
 end
 E_tr_m=mean(energy_tr_frame);
%% Total energy
energy_tot_frame=energy_rot_frame+energy_tr_frame;
E_tot_m=mean(energy_tot_frame);
Ti=E_tot_m/2.5/kB;