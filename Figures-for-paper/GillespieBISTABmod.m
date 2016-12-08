function [ mout,pout] = GillespieBISTABmod(N,par,totalreps,mstart, pstart,Output_Times);
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

P0=par(3);
MUM=par(2);
MUP = MUM;
n = par(1);
F = 10;

mout = zeros(totalreps,length(Output_Times));
pout = zeros(totalreps,length(Output_Times));

for reps=1:totalreps %parfor
% disp(reps)    
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
a1 = F*MUM*m;
a2 = F*MUM*p; 
a3 = F*N/(1 + ((p/N)/P0)^n);
a4 = F*N/(1 + ((m/N)/P0)^n);
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
                a1=F*(MUM)*(mn);
                a4 = F*N/(1 + ((mn/N)/P0)^n);             
            else
                if r2*a0<=a1
                    mn= m-1;
                    pn = p;
                    a1=F*MUM*(mn);
                    a4 = F*N/(1 + ((mn/N)/P0)^n);                 
                %elseif le(a1,r2*a0)*le(r2*a0,(a1+a2))
                elseif a1<=r2*a0 && r2*a0<=(a1+a2)
                    mn = m;
                    pn = p-1;
                    a2=F*MUP*pn;                    
                    a3 = F*N/(1 + ((pn/N)/P0)^n);
                elseif (a1+a2)<=r2*a0 && r2*a0<=(a1+a2+a3)
 
                    mn= m+1;
                    pn = p;
                    a1=F*MUM*(mn);
                    a4 = F*N/(1 + ((mn/N)/P0)^n);                     
                else

                    mn = m;
                    pn = p+1;                    
                    a2=F*MUP*pn;                    
                    a3 = F*N/(1 + ((pn/N)/P0)^n);                      
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