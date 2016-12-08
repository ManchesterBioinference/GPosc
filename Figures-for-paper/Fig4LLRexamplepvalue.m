load fig4

% fit models to data 



par1M = zeros(1,3);
par2M = zeros(1,4);
BICdiffM = zeros(1,1);

% load cell data for current experiment - loops through cells

i = 8;
y1 = dataNORMED(:,i);  

% fit OU and OUoscillatory models

[BICdiff, par1, par2] = getBICdiffRND(x,y1,Noise);
par1M(i,:) = par1;
par2M(i,:) = par2;
BICdiffM(i,:) = BICdiff;


%%

repeats = 2000;

[ synthOUhier1 ] = MakesynthOUHIERACHICAL(par1,repeats,x);

% [ BICdiffsynthTOT ] = BICdistDATAnew( synthOUhier1,time,par1TOT1,repeats,-4.5);
[ BICdiffM ] = BICdistDATAsynth( synthOUhier1,x,par1,repeats);
save('Fig4Nonoscexample')


%%


par1M = zeros(1,3);
par2M = zeros(1,4);
BICdiffM = zeros(1,1);

% load cell data for current experiment - loops through cells

i = 110;
y1 = dataNORMED(:,i);  

% fit OU and OUoscillatory models

[BICdiff, par1, par2] = getBICdiffRND(x,y1,Noise);
par1M(i,:) = par1;
par2M(i,:) = par2;
BICdiffM(i,:) = BICdiff;


%%

repeats = 2000;

[ synthOUhier1 ] = MakesynthOUHIERACHICAL(par1,repeats,x);

% [ BICdiffsynthTOT ] = BICdistDATAnew( synthOUhier1,time,par1TOT1,repeats,-4.5);
[ BICdiffM ] = BICdistDATAsynth( synthOUhier1,x,par1,repeats);
save('Fig4Oscexample')

%%

Thresh = BICdiff;
q1(i) = (sum(BICdiffsynthTOT>Thresh)/length(BICdiffsynthTOT));
