% N census from Hypergeom of biallelic haploid SNP alleles, full spectrum
% Based on NcHyper230810RandP2pass.m, spectrum of p each locus
% Simulates sampling inds&loci, uses jacknifed subsamples
% HypExp from HypExactExp230607.m, test (O-Ehyp) scaled 0-1 each locus
filename=input('outfile=','s');
%%%%%%%%%%%%%%%% Later make loops for samp etc
% Jacknife Sampling
samp=input('Largest fraction to sample, eg 0.1 to 0.9 = ');
jj=input('Prop of sample in each subsamp = ');

%  Increments Ptry, Ntry (hyps. for joint est. of Nc & p)
% 11 possible Ncensus values to trial (hyps for joint estimate Nc & p)
minNtry=input('minNtry to hypothesise '); maxNtry=3*minNtry;
Ninc=round((maxNtry-minNtry)/10); % fine increment, to give 11 Ntry values
Ntry=[minNtry:Ninc:maxNtry]; % array of hypothetical Nc values to trial
disp('N values to try'); Ntry
% Actual Nc value
Nit=input('Real Nit integer b/n min&maxN =');

% actual p value, called 'pp' for each locus, adjust to 'p' later
ppinc=0.05;% Increment from one pp value to next
ppcutter=[ppinc:ppinc:1-ppinc]; % the values
ppreps=10; % number of loci with each pp value
pp= sort(repmat(ppcutter,[1,ppreps])); % one value for each locus
L=length(pp); % Number of loci, index 'l'
for l=1:L
target=round((pp(l))*minNtry);p(l)=target/minNtry;%integral nbr targets
end

% hypothetical p values, Ptry; not same as p(l) above, esp. after jacknife
Ptry=[ppinc:ppinc:1-ppinc];  

%% Simulate Pop size Nit: individuals in population, each locus
preindit=zeros(L,Nit); % Rows Loci, Cols Inds
indit=zeros(L,Nit); % Rows Loci, Cols Inds
for l=1:L % popsetup-loci 1:L(row),inds (col - Nit  of these)
   targetit=round(p(l).*Nit); % nbr. target alleles in actual pop., locus 'l'
   preindit(l,1:targetit)=1; % set target inds to unity, rest still zero
   tempit=preindit(l,:);%vector of inds for locus 'l', length Nit
   idx = randperm(length(tempit));%idx:vec of shuffled indices of 'tempit'
   tempitperm = tempit(idx);%chooses inds (0 or 1) from 'tempit'
   indit(l,:)=tempitperm;%adds to matrix of L loci(row)* individs (col)
end % End locus loop to simulate population size Nit, L loci

% Loop for Ntry values, setting adaptive dimensions
%NB ensures n<Nit or n<Ntry, whichever is smaller
for Nindex=1:(length(Ntry)); % loop trial N's each locus, adaptive sampling
    %Set sample size n,jacknife subsampsize j(Nindex)
    if Nit<Ntry(Nindex); 
        n(Nindex)=round(Nit*samp); j(Nindex)=round(n(Nindex)*jj);
    else; 
        n(Nindex)=round(samp*Ntry(Nindex)); j(Nindex)=round(n(Nindex)*jj);
    end; % end adaptive sampling
    jmax(Nindex)=round(jj*round(maxNtry*samp)); % Common DIM all O & E vectors
end; % end Ntry loop for setting adaptive dimensions
    % Fixed DIMS, zereo fill to jmaxALL, for all "l", "Ntry".
    jmaxALL=max(jmax);
    HypExp=zeros(length(Ntry),length(Ptry),jmaxALL+1);
    propdetNit=zeros(L,jmaxALL+1);
    BIAS1=zeros(length(Ntry),length(Ptry),L);% BIAS for each locus  
% Loop for Ntry values, calc Hyp Exps for N,n,p
for Nindex=1:(length(Ntry)); % loop trial N's each locus,
    for pindex=1:(length(Ptry)); %loop Hyp Exp in sample of 'j', iTry230831n N loop
    parg=Ptry(pindex);Narg=Ntry(Nindex);narg=j(Nindex);%Arguments for HypExp next
    [HypExpTEMP1(1:(j(Nindex)+1))]=HypExactExp230607(parg,Narg,narg);
    HypExp(Nindex,pindex,(1:j(Nindex)+1))=HypExpTEMP1;%filled to jmax with zeros
    end; % end loop for Hyp Exp
end; % end Loop for Ntry values, calc Hyp Exps for N,n,p

% Main Loop for Ntry values
for Nindex=1:(length(Ntry)); % loop trial N's each locus,
    % Bin edges to tally targets-detected, increment 1; -0.5 to +(n+0.5)
    counter1=1;
    for e=-0.5:1:(j(Nindex)+0.5);%"j(Nindex)+1" bins, centres 0...j(Nindex)
    edges(counter1)=e;
    counter1=counter1+1;
    end; % end e loop 

    %% Major Locus loop begins
%% Jacknife sample actual population size Nit
for l=1:L; %sample loci Nit. Inds random in locus*popsize matrix
  %ploc(l)=(sum(indit(l,1:n)))/n;%allele prop locus "l", from sample of n
  for jndx=1:(n(Nindex)-j(Nindex)+1);% overlap jac of j(Nindex),loc'l'
  inditTEMP(l,(1:j(Nindex)))=indit(l,(jndx:(jndx+j(Nindex)-1)));   
  detNit(l,jndx)=sum(inditTEMP(l,1:j(Nindex)));%TrgtSampd row-loc,col-jack
  end; % end loop for 'n-j(Nindex)+1' subsamples each sample size j(Nindex)
end; % end loop for sampling over loci 'l'
 
% Tally results in 'freqdetNit', popsize Nit, & each locus
for l=1:L % loop loci for histogram data.
HistFreq=histcounts(detNit(l,:),edges);%in n-j(Nindex)+1 Jacknife reps
% row vec, bins number detected 0,1,...,j(Nindex)+1
PropNit=HistFreq./sum(HistFreq);
propdetNit(l,(1:j(Nindex)+1))=PropNit;%filled to jmaxALL with zeros
end; % end target-det; dims: locus, bins (0's,1's,...n detected)      
end; % end sampling fot that Ntry array

%% BIAS for each locus, & each hypothesised Ntry, Ptry
for l=1:L; % locus loop, each needs to examine all Ntrys
 for Nindex=1:(length(Ntry)); % loop trial N's each locus, 
  for pindex=1:(length(Ptry)); % begin joint estimate p & Nc, loop p values
    % AbsBias, Absolute value:Obs-HypergeoExps, n-1 inds  (count 0,1,...n)
    HypExpTEMP2=(zeros(1,jmaxALL+1))';
    aa=HypExp(Nindex,pindex,(1:jmaxALL+1));
    HypExpTEMP2=(squeeze(aa))'; % row vec length j(Nindex)
    PropTEMP=propdetNit;
    BiasTEMP1=abs(PropTEMP-HypExpTEMP2);%fill to jmax with 0
    BiasTEMP2=sum(BiasTEMP1(l,:)); % abs bias over all bins, this locus
    BIAS1(Nindex,pindex,l)=BiasTEMP2;% AbsBias all bins, this locus 
  end; % end pindex loop
 end; % End loop for simulating and testing Ntry values.
end; %end locus loop.

% Joint Pest-Nest diagnosis each locus, then Nest averaged over all loci
%Diagnose p,N this locus, search min BIAS1, index gives Ptry Ntry indices
for l=1:L;
locBIAS=squeeze(BIAS1(:,:,l)); %bias for this locus (3rd dim)
Min=min(locBIAS,[],"all");
[row,col]=find(locBIAS==Min); % shows row and col of min value
Nest(l)=mean(Ntry(row)); Pest(l)=mean(Ptry(col));% Av if>1 val
end; % end locus loop
% Over loci, AVE and SE of Nest
AveNest=mean(Nest); SeNest=(std(Nest))/sqrt(L);

% For investigating loci where Pest Biased up so Nest Biased down.
NestBias(:)=(Nest(:)-Nit)./Nit;
PestBias(:)=(Pest(:)-p(:))./p(:);
SDdetNit=(std(detNit,1,2))';
AvedetNit=(mean(detNit,2))';
CVdetNit=((SDdetNit(:))./(AvedetNit(:)))';
% InvCVdetNit=((AvedetNit(:))./(SDdetNit(:)))';

%% Now do Nest diagnosis with low CV loci only
LocNbr(1,1:L)=[1:L];% serialnumbers for loci
CVNest=([CVdetNit(1:L);Nest(1:L);LocNbr(1:L)])';%CVcol1,NestCol2;LocNbrCol3
CVNestSort= (sortrows(CVNest,1))';%sorts by CV, all rows (CV,Nest,LocNbr)
NestLoCV=CVNestSort(2,(1:round(L/10)));%Nest row,lowest 10% of CVs    
CVLoLoci=CVNestSort(3,(1:round(L/10)));% locus numbers used in AveNestLoCV
AveNestLoCV=mean(NestLoCV);
SDNestLoCV=std(NestLoCV);

% Plot NestBias PestBias vs CV,skew, etc of detNit(l,jndx)
subplot(2,2,1); % NestBias vs CVdetNit
xlabel('CVdetNit');
ylabel('NestBias');
p1=plot(CVdetNit,NestBias);

subplot(2,2,2); % PestBias vs CVdetNit 
xlabel('CVdetNit');
ylabel('PestBias');
p1=plot(CVdetNit,PestBias);

%subplot(2,2,3); % NestBias vs InvCVdetNit
%xlabel('InvCVdetNit');
%ylabel('NestBias');
%p1=plot(InvCVdetNit,NestBias);

%subplot(2,2,4); % PestBias vs InvCVdetNit 
%xlabel('InvCVdetNit');
%ylabel('PestBias');
%p1=plot(InvCVdetNit,PestBias);
%% PLOTTING
% SumBIAS2=squeeze(sum(BIAS2,3)); % Sum bias over all loci (3rd dim)
% MeanBIAS2=SumBIAS2./L; % Mean bias over all loci, for plotting
% Y=Ntry2; X=Ptry2; % Axes for contour plot
% contour(X,Y,MeanBIAS2,'showtext','on');

save(filename,'p','L','Nit','Ntry','maxNtry','Ptry','samp','jj','n','j','indit','propdetNit','HypExp','BIAS1','Nest','Pest','AveNest','SeNest','CVNest','CVLoLoci','AveNestLoCV','SDNestLoCV','NestLoCV','CVLoLoci','AveNestLoCV','SDNestLoCV')
clearvars -except filename

%% PASS2:  Uses the alternative allele
load(filename,'p','L','Nit','Ntry','maxNtry','Ptry','samp','jj','indit')
inditTEMP=indit;
indit=(inditTEMP-1).*(-1);% swaps zeros and ones 

% Loop for Ntry values, setting adaptive dimensions
%NB ensures n<Nit or n<Ntry, whichever is smaller
for Nindex=1:(length(Ntry)); % loop trial N's each FPestlocus, adaptive sampling
    %Set sample size n,jacknife subsampsize j(Nindex)
    if Nit<Ntry(Nindex); 
        n(Nindex)=round(Nit*samp); j(Nindex)=round(n(Nindex)*jj);
    else; 
        n(Nindex)=round(samp*Ntry(Nindex)); j(Nindex)=round(n(Nindex)*jj);
    end; % end adaptive sampling
    jmax(Nindex)=round(jj*round(maxNtry*samp)); % Common DIM all O & E vectors
end; % end Ntry loop for setting adaptive dimensions
    % Fixed DIMS, zereo fill to jmaxALL, for all "l", "Ntry".
    jmaxALL=max(jmax);
    HypExp=zeros(length(Ntry),length(Ptry),jmaxALL+1);
    propdetNit=zeros(L,jmaxALL+1);
    BIAS1=zeros(length(Ntry),length(Ptry),L);% BIAS for each locus  
% Loop for Ntry values, calc Hyp Exps for N,n,p
for Nindex=1:(length(Ntry)); % loop trial N's each locus,
    for pindex=1:(length(Ptry)); %loop Hyp Exp in sample of 'j', iTry230831n N loop
    parg=Ptry(pindex);Narg=Ntry(Nindex);narg=j(Nindex);%Arguments for HypExp next
    [HypExpTEMP1(1:(j(Nindex)+1))]=HypExactExp230607(parg,Narg,narg);
    HypExp(Nindex,pindex,(1:j(Nindex)+1))=HypExpTEMP1;%fill jmaxALL+1 with0
    end; % end loop for Hyp Exp
end; % end Loop for Ntry values, calc Hyp Exps for N,n,p

% Main Loop for Ntry values
for Nindex=1:(length(Ntry)); % loop trial N's each locus,
    % Bin edges to tally targets-detected, increment 1; -0.5 to +(n+0.5)
    counter1=1;
    for e=-0.5:1:(j(Nindex)+0.5);%"j(Nindex)+1" bins, centres 0...j(Nindex)
    edges(counter1)=e;
    counter1=counter1+1;
    end; % end e loop 

    %% Major Locus loop begins
%% Jacknife sample actual population size Nit
for l=1:L; %sample loci Nit. Inds random in locus*popsize matrix
  %ploc(l)=(sum(indit(l,1:n)))/n;%allele prop locus "l", from sample of n
  for jndx=1:(n(Nindex)-j(Nindex)+1);% overlap jac of j(Nindex),loc'l'
  inditTEMP(l,(1:j(Nindex)))=indit(l,(jndx:(jndx+j(Nindex)-1)));   
  detNit(l,jndx)=sum(inditTEMP(l,1:j(Nindex)));%TrgtSampd row-loc,col-jack
  end; % end loop for 'n-j(Nindex)+1' subsamples each sample size j(Nindex)
end; % end loop for sampling over loci 'l'
 
% Tally results in 'freqdetNit', popsize Nit, & each locus
for l=1:L % loop loci for histogram data.
HistFreq=histcounts(detNit(l,:),edges);%in n-j(Nindex)+1 Jacknife reps
% row vec, bins number detected 0,1,...,j(Nindex)+1
PropNit=HistFreq./sum(HistFreq);
propdetNit(l,(1:j(Nindex)+1))=PropNit;%filled to jmaxALL+1 with zeros
end; % end target-det; dims: locus, bins (0's,1's,...n detected)      
end; % end sampling fot that Ntry array

%% BIAS for each locus, & each hypothesised Ntry, Ptry
for l=1:L; % locus loop, each needs to examine all Ntrys
 for Nindex=1:(length(Ntry)); % loop trial N's each locus, 
  for pindex=1:(length(Ptry)); % begin joint estimate p & Nc, loop p values
    % AbsBias, Absolute value:Obs-HypergeoExps, n-1 inds  (count 0,1,...n)
    HypExpTEMP2=(zeros(1,jmaxALL))';
    aa=HypExp(Nindex,pindex,(1:jmaxALL+1));
    HypExpTEMP2=(squeeze(aa))'; % row vec length j(Nindex)
    PropTEMP=propdetNit;
    BiasTEMP1=abs(PropTEMP-HypExpTEMP2);%fill to jmax with 0
    BiasTEMP2=sum(BiasTEMP1(l,:)); % abs bias over all bins, this locus
    BIAS1(Nindex,pindex,l)=BiasTEMP2;% AbsBias all bins, this locus 
  end; % end pindex loop
 end; % End loop for simulating and testing Ntry values.
end; %end locus loop.

% Joint Pest-Nest diagnosis each locus, then Nest averaged over all loci
%Diagnose p,N this locus, search min BIAS1, index gives Ptry Ntry indices
for l=1:L;
locBIAS=squeeze(BIAS1(:,:,l)); %bias for this locus (3rd dim)
Min=min(locBIAS,[],"all");
[row,col]=find(locBIAS==Min); % shows row and col of min value
Nest(l)=mean(Ntry(row)); Pest(l)=mean(Ptry(col));% Av if>1 val
end; % end locus loop
% Over loci, AVE and SE of Nest
AveNest=mean(Nest); SeNest=(std(Nest))/sqrt(L);

% For investigating loci where Pest Biased up so Nest Biased down.
NestBias(:)=(Nest(:)-Nit)./Nit;
PestBias(:)=(Pest(:)-p(:))./p(:);
SDdetNit=(std(detNit,1,2))';
AvedetNit=(mean(detNit,2))';
CVdetNit=((SDdetNit(:))./(AvedetNit(:)))';
% InvCVdetNit=((AvedetNit(:))./(SDdetNit(:)))';

%% Now do Nest diagnosis with low CV loci only
LocNbr(1,1:L)=[1:L];% serialnumbers for loci
CVNest=([CVdetNit(1:L);Nest(1:L);LocNbr(1:L)])';%CVcol1,NestCol2;LocNbrCol3
CVNestSort= (sortrows(CVNest,1))';%sorts by CV, all rows (CV,Nest,LocNbr)
NestLoCV=CVNestSort(2,(1:round(L/10)));%Nest row,lowest 10% of CVs    
CVLoLoci=CVNestSort(3,(1:round(L/10)));% locus numbers used in AveNestLoCV
AveNestLoCV=mean(NestLoCV);
SDNestLoCV=std(NestLoCV);

%% End PASS2 with allele 2

%% change names to avoid overwrite of PASS1
p2=p;L2=L;Nit2=Nit;Ntry2=Ntry;maxNtry2=maxNtry;Ptry2=Ptry;
samp2=samp;jj2=jj;n2=n;j2=j;indit2=indit;
propdetNit2=propdetNit;HypExp2=HypExp;BIAS1o2=BIAS1;
Nest2=Nest;Pest2=Pest;AveNest2=AveNest;SeNest2=SeNest;CVNest2=CVNest;
CVLoLoci2=CVLoLoci;AveNestLoCV2=AveNestLoCV;SDNestLoCV2=SDNestLoCV;
NestLoCV2=NestLoCV;

%% Now make NestBoth by combined av Nest 1 +2
NestLoCV2=NestLoCV;
NestLoCVTEMP=load(filename,'NestLoCV');%gets back Nest for other allele (overwritten)
workaround=struct2cell(NestLoCVTEMP);NestLoCV1=[workaround{:}];%workaround
NestBothLoCV=[NestLoCV1,NestLoCV2];% concatenates NestLoCV a1 & a2
AveNestBothLoCV=mean(NestBothLoCV);
SDNestBothLoCV=std(NestBothLoCV);

%% Binomial better fit than Hyp? If so, need bigger sample genotyped
% for most likely values AveNestLoCV and AvePest(1:L)=(Pest+Pest2)/2

load(filename,'Pest');% get back Pest from PASS1
Pest2reflect=1-Pest2; % Put Pest2 in same direction as Pest2
AvePest(:)=(Pest(:)+Pest2reflect(:))./2;% Av p estimate for each locus
nNest=round(AveNestBothLoCV*samp);jNest=round(nNest*jj);%appropriate n,j Value
BinomExp=zeros(L,jNest+1);%fill to jNest+1 with zeros
HypExpFinal=zeros(L,jNest+1);%fill to jNest+1 with zeros
HypExpFinalTEMP=zeros(L,jNest+1);%fill to jNest+1 with zeros

for l=1:L
  % Binomial Exps for subsample of 'jNest',AvePest(l)
  bins=[0:1:jNest];% bins for BinomExps
  parg=AvePest(l);narg=jNest;% Arguments for next line
  [BinomExpTEMP1(1:jNest+1)]=binopdf(bins,narg,parg);
  BinomExp(l,(1:jNest+1))=BinomExpTEMP1;%1*jmaxALL, trailing zeros
  % Hyp Exps for AveNestLoCV,jNest,AvePest(l)
  parg=AvePest(l);narg=jNest;Narg=round(AveNestLoCV);%Arguments for next
[HypExpFinalTEMP(l,(1:(jNest+1)))]=HypExactExp230607(parg,Narg,narg);%1*jNest+1
end; % end locus loop
% Fix numbers so low they are NaN ("not a number")
TF = isnan(HypExpFinalTEMP);% "1" wherever value so low it is NaN, others"0"
HypExpFinal=HypExpFinalTEMP;
HypExpFinal(TF)=0;

% Get final BIAS
load(filename,'propdetNit');%propdetNit is L*jmaxALL
propdetTRIM=propdetNit(:,1:jNest+1); % trim trailing zeros
BIASBinom=abs(propdetTRIM-BinomExp); %DIM row 1:L, Col 1:jNest+1
BIASbinVec=reshape(BIASBinom,[1,numel(BIASBinom)]);%for t-test
BIASbin=sum(BIASBinom,"all");% NB lower case 'b' sum all bins & all loci
BIASHyp=abs(propdetTRIM-HypExpFinal);
BIAShypVec=reshape(BIASHyp,[1,numel(BIASHyp)]);%linear for t-test
BIAShyp=sum(BIASHyp,"all");% NB lower case 'h' sum all bins & all loci

% Useless test whether Hypergeom fits jacknife data better than binomial
[h,psig] = ttest(BIAShypVec,BIASbinVec);
SamSize=BIAShyp-BIASbin; % positive if hypergeom fits better than binomial 
%if psig<=0.05
%    if SamSize>0
%  disp('Need bigger sample, SamSize= '); SamSize
%  elseif SamSize<0 disp ('Hypergeom fits better than Binomial, sample OK');
% end; % end if SamSize
% end; % end if psig

save(filename,'p2','L2','Nit2','Ntry2','maxNtry2','Ptry2','samp2','jj2','n2','j2','indit2','propdetNit2','HypExp2','BIAS1o2','Nest2','Pest2','AveNest2','SeNest2','CVNest2','CVLoLoci2','AveNestLoCV2','SDNestLoCV2','NestLoCV2','CVLoLoci2','AveNestLoCV2','SDNestLoCV2','NestBothLoCV','AveNestBothLoCV','SDNestBothLoCV','BIASBinom','BIASHyp','psig','SamSize','NestBothLoCV','AveNestBothLoCV','SDNestBothLoCV',"-append")

%% Full Binomial if needed
% Binomial Exps for subsample of 'j(Nindex)'
%BinomExp=zeros(length(Ntry),length(Ptry),jmaxALL);
%for Nindex=1:(length(Ntry)); %loop trial N's
%for pindex=1:(length(Ptry)); %loop trial P's
%    bins=[0:1:j(Nindex)];% bins for BinomExps
%    parg=Ptry(pindex);narg=j(Nindex);% Arguments for next line
%    [BinomExpTEMP1(1:(j(Nindex)+1))]=binopdf(bins,narg,parg);
%    BinomExp(Nindex,pindex,(1:j(Nindex)+1))=BinomExpTEMP1;%filled to jmax with zeros
%end; % end loop for BinomExp Ptry values
%end; % end loop for BinomExp Ntry values