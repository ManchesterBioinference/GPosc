function [ mout,pout] = GillespieTiming_Mod9nodelay(N,par,totalreps,mstart, pstart,Output_Times);
%DELAYGILLESPIE Runs Gillespie algorithm with delay processes. Returns a
% value of time (dist)
%unpack parms
% par(1) = 10;
% par(2) = 4.1;
% par(3) = 0.03;
% par(4) = 0.03;
% par(5) = 1;
% par(6) = 1;
% par(7) = 18.7;

P0=par(1);
NP=par(2);
MUM = par(3);
MUP= par(4);
ALPHAP=par(6);
ALPHAM=par(5);
tau = par(7);

mout = zeros(totalreps,length(Output_Times));
pout = zeros(totalreps,length(Output_Times));

for reps=1:totalreps %parfor
disp(reps)    
j_t_next = 1;
    
iter = 1;
rlist = [];

% sets initial values
m=mstart;
mn=0;
p=pstart;
pn=0;
t=0;tn=0;
rlist=[];

%When rlist is empty, there is a different process:
%We follow this process until rlist is no longer empty
a1 = MUM*m;
a2 = MUP*p; 
a3 = ALPHAP*m;
a4 = N*ALPHAM/(1 + ((p/N)/P0)^NP);
while t<Output_Times(end) % runs until p hits a threshold            
            a0=a1+a2+a3+a4;
            r1=rand(1);
            r2=rand(1);
            dt=(1/a0)*log(1/r1);
            if numel(rlist)>0 && t<=rlist(1) && rlist(1)<=(t+dt)
            %if le(t(iter),rlist(1))*le(rlist(1),(t(iter)+dt))
                mn= m+1;
                pn = p;
                tn = rlist(1);
                rlist(1)=[];
                a1=(MUM)*(mn);
                a3=ALPHAP*mn;                
            else
                if r2*a0<=a1
                    mn= m-1;
                    pn = p;
                    a1=MUM*(mn);
                    a3=ALPHAP*mn;                    
                %elseif le(a1,r2*a0)*le(r2*a0,(a1+a2))
                elseif a1<=r2*a0 && r2*a0<=(a1+a2)
                    mn = m;
                    pn = p-1;
                    a2=MUP*pn;                    
                    a4=N*ALPHAM/(1 + ((pn/N)/P0)^NP);
                elseif (a1+a2)<=r2*a0 && r2*a0<=(a1+a2+a3)
                    mn = m;
                    pn = p+1;                    
                    a2=MUP*pn;                    
                    a4=  N*ALPHAM/(1 + ((pn/N)/P0)^NP);                
                else
                    mn= m+1;
                    pn = p;
                    a1=MUM*(mn);
                    a3=ALPHAP*mn;                    
                end
                tn = t+dt;
            end
        iter = iter + 1;
        t=tn;
        m =mn;
        p =pn;        
        while j_t_next<=length(Output_Times)&&t>Output_Times(j_t_next)
            mout(reps,j_t_next) = m;
            pout(reps,j_t_next) = p;            
            j_t_next=j_t_next+1;
        end
end 
   
end                   
                                     
end

