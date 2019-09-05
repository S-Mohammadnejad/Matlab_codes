%Extract the final velocities in the H2 box and verify if they follow Maxwell distribution
clear
close all
clc
%% Input parameters
%% Properties
kB=1.38064852*10^(-23); %[SI unit]
m_Ar=39.948/1000/(6.02*10^23); %molecular mass [kg]
m_H=4.002/1000/(6.02*10^23); %molecular mass [kg]
m_H2=2*1.0079/1000/(6.02*10^23); %molecular mass [kg]
mg=m_H2; %<<decide which TAC to be calculated: argon or hydrogen
mu=m_H2/4; %reduced mass for H2 molecule
%<<change temperatures
Tw=298; %surface temperature [K]
Tg=295; %gas molecule shot/initial temperature [K]
r_cut=10; % cutt-off radius in ReaxFF
v_2_w=4*kB*Tw/mg;
v_mp_Tw=sqrt(2*kB*Tw/mg);
E_w_tr=2*kB*Tw;
E_w_ro=kB*Tw;
E_w_tot=3*kB*Tw;
E_g_tr=1.5*kB*Tg;
E_g_ro=kB*Tg;
E_g_tot=2.5*kB*Tg;
sigma=sqrt(kB*Tg/mg);
%% Read files
logfiles=dir('*.logfile');
num_logs=length(logfiles);
rxxmols=dir('*.rxxmol');
num_rxxmols=length(rxxmols);
%cd(oldFolder);
[ num_frame, Models, num_try ] = Read_sequence( logfiles,num_logs,rxxmols,num_rxxmols );
fname=logfiles(1).name;
[time_span,output_interval]=logfile_scan_H2_box(fname);
v_ref=0:0.01:10000;
f_v=4*pi*(mg/(2*pi*kB*Tg))^1.5.*v_ref.^2.*exp(-0.5*mg*v_ref.^2./(kB*Tg)); %Maxwell distribution
f_v2=0.5*(mg/(kB*Tg))^2*v_ref.^3.*exp(-mg.*v_ref.^2./(2*kB*Tg)); %Incident velocity distribution

m=zeros(1,num_try);%%%%%%in eaxh rxxmol we want to know how many atoms are
for i=1:num_try
    m(i)=Models{i}(num_frame(i)).MolNo;
end
%%
bar=20; %number of bars
%%%%%%%%%%%%%%%%%%%%
%%%-------WE COMPUTE THE C.O.M VELOCITY FOR THE LAST FRAME OF TRIJECTORY
%%%%%%%%%%%%%%%%%%%%%
vc=zeros(3,1);
a=1;
c=0;
for i=1:num_try
    for j=a:2:(i*m(i))
        vc(:,j)=(Models{i}(num_frame(i)).Velocity(j-c,:)+Models{i}(num_frame(i)).Velocity(j+1-c,:))/2; %he is computing mean just for the last frame
        if j==i*m(i)-1
            a=i*m(i)+1;
            c=i*m(i);
            break
        end
    end
end
vc=vc(:,any(vc)); %deleting zero elements in the matrix
abs_vc=sqrt(dot(vc,vc)); %dot computes the dot product alonge columns
%%
figure(1)
histogram(abs_vc,bar,'Normalization','pdf')
hold on
plot(v_ref,f_v);
hold off
xlabel('v_{c} (m/s)');
ylabel('f(v_{c})')
legend('Sampled velocities','Maxwell distribution','Location','Best');
%%
%------------converting velocities---------------------
%criterion based on velocity at center of mass
H2_vel=zeros(3,1);
b=1;
c=0;
for j=1:num_try
for i=b:2:(j*m(j))
    vc_z=(Models{j}(num_frame(j)).Velocity(i-c,3)+Models{j}(num_frame(j)).Velocity(i+1-c,3))/2; %we use last frame as the reference
    %normal component of velocity at center of mass [m/s]
    if vc_z<0
        H2_vel(:,i)=Models{j}(num_frame(j)).Velocity(i-c,:);
        H2_vel(:,i+1)=Models{j}(num_frame(j)).Velocity(i+1-c,:);
    end
    if vc_z>0
        H2_vel(1,i)=Models{j}(num_frame(j)).Velocity(i-c,1);
        H2_vel(2,i)=Models{j}(num_frame(j)).Velocity(i-c,2);
        H2_vel(3,i)=-Models{j}(num_frame(j)).Velocity(i-c,3);
        H2_vel(1,i+1)=Models{j}(num_frame(j)).Velocity(i+1-c,1);
        H2_vel(2,i+1)=Models{j}(num_frame(j)).Velocity(i+1-c,2);
        H2_vel(3,i+1)=-Models{j}(num_frame(j)).Velocity(i+1-c,3);
    end
    if vc_z==0
        continue
    end
    if i==j*m(j)-1
        b=j*m(j)+1;
        c=j*m(j);
        break
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-----I add this

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=find(H2_vel(3,:));
H2_vel=H2_vel(:,n);
%%
%-------sampling based on crossing the plane---------------
%-
%-----it is nice method for computing COM velocity
vc_z=zeros(1,length(H2_vel(3,:)));
for j=1:2:length(H2_vel(3,:))
    vc_z(j)=(H2_vel(3,j)+H2_vel(3,j+1))/2;
end
vc_z=vc_z(vc_z~=0); %normal velocities at center of mass its 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vc_z_A=zeros(1,length(vc_z));
H2_vel_A=zeros(3,2*length(vc_z));
H2_vel_filter=zeros(1,2*length(vc_z));
dt=0.25; %random time step, trivial [s]
dz_max=max(-vc_z)*dt; %maximum distance
dz=-vc_z.*dt;

for k=1:length(vc_z)
    R=rand(1);
    
    if R*dz_max<dz(k)
        n=2*k-1;
        
        vc_z_A(k)=vc_z(k);
        
        H2_vel_A(:,n)=H2_vel(:,n);      
        H2_vel_A(:,n+1)=H2_vel(:,n+1);
        H2_vel_filter(n)=n;
        H2_vel_filter(n+1)=n+1;
    end
end
nc_z=find(H2_vel_A(1,:)==0);
vc_z_A=vc_z_A(vc_z_A~=0).*10^10; %[Angstrom/s]\

H2_vel_less=zeros(3,length(nc_z)); %these are velocity components which do not satisfy above condition
for i=1:length(nc_z)
    H2_vel_less(:,i)=H2_vel(:,nc_z(i)).*10^10; %[Angstrom/s];
end
%%%%%%
%%%%%
H2_vel_A=H2_vel_A(:,any(H2_vel_A)).*10^10; %[Angstrom/s]
H2_vel_final=[H2_vel_A H2_vel_less];
%% save the velocities according to ADF format requirements
%
fid=fopen('H2_vel_295.txt','w+');
fprintf(fid,'v[%d]="%24.17E%24.16E%24.16E"\r\n',[1:length(H2_vel_final(3,:));H2_vel_final]);
fclose(fid);
%%
save('H2_vel.mat','H2_vel','-v7.3')
save('H2_vel_A.mat','H2_vel_A','-v7.3')
save('H2_vel_less.mat','H2_vel_A','-v7.3')
save('variables.mat','m','num_frame','num_try','output_interval')
%% Plot
figure(5)
v_z=-5000:1:0;
f_vz=-v_z.*exp(-v_z.^2./(2*sigma^2))./(sigma^2); %theoretical distribution for normal component (vz)
%histogram(vc_z_A./10^10,bar,'Normalization','pdf') %actual distribution (vz)
histogram(vc_z,20,'Normalization','pdf')
hold on
plot(v_z,f_vz)
hold off
xlabel('v_{cz} (m/s)');
ylabel('f(v_{cz})')
legend('Sampled velocities','Rayleigh distribution','Location','Best');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
v_z=-5000:1:0;
f_vz=-v_z.*exp(-v_z.^2./(2*sigma^2))./(sigma^2); %theoretical distribution for normal component (vz)
histogram(vc_z_A./10^10,bar,'Normalization','pdf') %actual distribution (vz)
hold on
plot(v_z,f_vz)
hold off
xlabel('v_{cz} (m/s)');
ylabel('f(v_{cz})')
legend('Sampled velocities','Rayleigh distribution','Location','Best');
figure(3)
vc_x_A=zeros(1,length(H2_vel_A(1,:)));
for j=1:2:length(H2_vel_A(1,:))
    vc_x_A(j)=(H2_vel_A(1,j)+H2_vel_A(1,j+1))/2;
end
vc_x_A=vc_x_A(vc_x_A~=0);
v_x=-5000:1:5000;
f_vx=(mg/(2*pi*kB*Tg))^0.5.*exp(-mg*v_x.^2./(2*kB*Tg)); %theoretical distribution for tangential component (vx)
histogram(vc_x_A./10^10,bar,'Normalization','pdf') %actual distribution (vx)
hold on
plot(v_x,f_vx)
hold off
xlabel('v_{cx} (m/s)');
ylabel('f(v_{cx})')
legend('Sampled velocities','Maxwell distribution','Location','Best');
figure(4)
vc_A=zeros(3,length(H2_vel_A(1,:)));
for j=1:2:length(H2_vel_A(1,:))
    vc_A(:,j)=(H2_vel_A(:,j)+H2_vel_A(:,j+1))/2/10^10;
end
vc_A=vc_A(:,any(vc_A));
abs_vc_A=sqrt(dot(vc_A,vc_A));
histogram(abs_vc_A,bar,'Normalization','pdf')
hold on
plot(v_ref,f_v,v_ref,f_v2);
hold off
xlabel('v_{c} (m/s)');
ylabel('f(v_{c})')
legend('Sampled velocities','Maxwell distribution','Incident velocity distribution','Location','Best');








