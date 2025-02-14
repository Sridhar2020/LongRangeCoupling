%%%LinkBoundary.m%%%

%%%Author S.Sridhar (dharmails@gmail.com)
%%%Affl: University of Sheffield, UK

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%This program generates myocyte-fibroblast links 
%%%that couple fibroblasts with both local and distal
%%%myocytes. The fibroblasts are locally attached
%%%to a region of myocytes in a circular scar and distally
%%%attached to a region of myocytes around the scar.
%%%The input parameters include 
%%% a) lam - the Poisson parameter determing the number of links
%%%    attached to a given fibroblast unit
%%% b) np  - the number of fibroblast units attached to the scar
%%% c) FDist,FDist1,FDist2 - the Euclidean distance between myocytes
%%%    upto which fibroblasts can couple myocytes.
%%% d) N*N - dimensions of the myocyte grid
%%% e) rad1 and rad2 - radius of scar and border zone
%%%    used to generate random M-F links
%%%   
    
%%% The output is a sparse adjacency matrix of size (N*N, np) with entries 1 or 0.
%%% The entries with 1 correspond to links while 0 implies no links
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

tic
%Set input parameters
np = 30000;lam=50;N=400;
FDist=10;FDist1=5;FDist2=0; %Euclidean distance in terms of grid points
rad1 = 100;rad2=80;
cent = 200; %N/2 
%Creating random links - Poisson distribution
f2=poissrnd(lam,N*N,1);%each fibroblast has a random no. of connections 
f2= f2(f2~=0); %get the fibroblasts with non-zero links 
%Fix np sites 
f1=f2(90001:120000); %in this case 30000 sites as np=30000
LL1 =[]; x1=[]; y1=[]; LLL1=[];
%Initialise adjacency matrices
adj = sparse(N*N,np);adj1=sparse(N*N,np);adj2=sparse(N*N,np);

%Identiying region of scar and border zone
for i = 1:N 
    for j = 1:N
        
        if ((sqrt((i-cent).^2 + (j-cent).^2) -rad1 )< 0)        
           fibbord((i-1)*N+j)=1;
        else
            fibbord((i-1)*N+j)=0;
        end
        
        if ((sqrt((i-cent).^2 + (j-cent).^2) -rad2) < 0)        
            fibscar((i-1)*N+j)=1;
        else
            fibscar((i-1)*N+j)=0;
        end
    end
end
pop1=find(fibbord); %myocytes in border zone with links to fibrobast
pop2=find(fibscar); %myocytes in scar with links to fibroblast

%Determine the adjacency matrix entries for each of the fibroblast
 for k = 1:np 
     
     %uniform random sample of which myocytes to be linked to the fibroblast
     x = randsample(pop1,f1(k)); % myocytes to be linked to fibroblast in border zone
     y=randsample(pop2,f1(k)); %myocytes to be linked to fibroblast in scar
     x1=[x1 x];y1=[y1 y];
    
     %If the fibroblast has only one link (from the Poisson distribution)
     %it is connected to the myocyte in the scar (chosen from randsample)
     if(length(x)==1)
       adj(x,k) = 0;
       adj1(x,k)=0;
       adj2(x,k)=0;
     end
     if(length(y)==1)
       adj(y,k) = 1;
       adj1(y,k)=1;
       adj2(y,k)=1;
     end

    %If the fibroblast has more than one link, it is connected to
    %proximal myocyte on grid and distal myocyte provided
    %condition on distance (FDist) is satisfied.
      if(length(x) > 1)
       adj(x(1),k)=0;adj1(x(1),k)=0;adj2(x(1),k)=0;
       for i = 1:length(x)
         j1(i)=mod(x(i),N);i1(i)=1+(x(i)-j1(i))/N;
       end
       for i = 1:length(x)-1 
        LL(i)=sqrt((i1(1)-i1(i+1)).^2 + (j1(1) - j1(i+1)).^2);
        LL1= [LL LL1];
         if (LL(i) > FDist || LL(i) == 0)
             adj(x(i+1),k)=0;
         elseif(LL(i) <= FDist && LL(i) > 0)
             adj(x(i+1),k)=1;
         end
         if (LL(i) > FDist1 || LL(i) == 0)
             adj1(x(i+1),k)=0;
         elseif(LL(i) <= FDist1 && LL(i) > 0)
             adj1(x(i+1),k)=1;
         end
         if (LL(i) > FDist2 || LL(i) == 0)
             adj2(x(i+1),k)=0;
         elseif(LL(i) <= FDist2 && LL(i) > 0)
             adj2(x(i+1),k)=1;
         end

       end
     end
     if(length(y) > 1)
       adj(y(1),k)=1;adj1(y(1),k)=1;adj2(y(1),k)=1; 
       for i = 1:length(y)-1 
             adj(y(i+1),k)=0;
             adj1(y(i+1),k)=0;
             adj2(y(i+1),k)=0;
       end
    end
 end
%Save output adjacency matrices for all 3 maximum distance of M-F links
 save('LinkMatrix_Size400_ScarRad80_and_BorderRad100_Lam50_Np30000.mat','adj','adj1','adj2','np','f1','FDist','FDist1','FDist2');
 toc

