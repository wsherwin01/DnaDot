% N census from Hypergeom of diploid genotypes
% Based on NcHyper240129compilable

%% Choice Inputs etc
% Jacknife Sampling: jj is proportion of sample in each subsample eg 0.6
load Choices.txt; % Cell A1 is jj, Cell A2 is minNtry
jj=Choices(1,1);
%  Increments Ptry, Ntry (hypotheses for joint est. of Nc & p)
% 11 possible Ncensus values to trial (hyps for joint estimate Nc & p)
minNtry=2*Choices(1,2); maxNtry=3*minNtry;%diploid data
Ninc=round((maxNtry-minNtry)/10); % fine increment, to give 11 Ntry values
Ntry=[minNtry:Ninc:maxNtry]; % array of hypothetical Nc values to trial
% hypothetical p values to trial
ppinc=0.05;% Increment from one pp value to next
Ptry=[ppinc:ppinc:1-ppinc];

%% Data Input
load Genotypes.txt;%row=individual, column-pair=locus
Input2=Genotypes'; % transposed col=individual, row-pair=locus
RowCol=size(Input2); TwiceLoc=RowCol(1); Inds=RowCol(2);
% transpose and concatenate each locus: row locus, col ind
Input3=Input2;Input4=Input2; %temps to remove every second row from
Input3(1:2:TwiceLoc,:)=[];
Input4(2:2:TwiceLoc,:)=[];
preindit=[Input3,Input4];
%preindit2=preindit; preindit2(preindit2==0)=1000;% zero length coded '1000'
% gives preindit row=loci, col=inds*2 for alleles 1&2
locind=size(preindit);
L1=locind(1); n=(locind(2))/2;% nbr true loci 'L1', total sample size 'n'
indit=preindit; % assumes already randomised by user

%% convert loci to ""loci"" - one for each allele at each locus
L=0; % Cumulate number ""loci"" =Sum(allele-nbr) over loci, in these loops
for l1=1:L1; % loop over true loci, get 'L' = nbr of ""loci""
   alcode = unique(indit(l1,:));%allele sizes for locus 'l1'
   alcodelen(l1)=length(alcode); % how many allele sizes
   for al=1:alcodelen(l1); %loop over alleles within this true locus
      LTEMP=L+al; % treats each allele as separate ""locus"
      indal(LTEMP,:)=indit(l1,:)==alcode(al);%"locus"",target coded '1'
   end;% end loop over alleles within this true locus
   L=LTEMP; %treats +/- each allele as separate ""locus""
end; % end true locus loop; indal row ""locus""(L),col ""individual"" (2*n)
indalTEMP1=indal; indalq=(indalTEMP1-1).*(-1);% swaps zeros and ones
indalfold=[indal;indalq]; % contains target and (1-target) alleles

%  allele proportions for each ""locus"", incl target & (1-target) 'p(l)'
for l=1:L; % ""locus"" loop
    sumTEMP=sum(indalfold(l,:));nTEMP=2*n;
    p(l)=sumTEMP/nTEMP;
end; % end ""locus"" loop

j=round(2*n*jj);
% ended data input and transform

% Fixed DIMS, zero fill to jmaxALL, for all "l", "Ntry".
    HypExp=zeros(length(Ntry),length(Ptry),j+1);
    propdetNit=zeros(L,j+1);
    BIAS1=zeros(length(Ntry),length(Ptry),L);% BIAS for each locus
% Loop for Ntry values, calc Hyp Exps for N,n,p
for Nindex=1:(length(Ntry)); % loop trial N's each locus,
    for pindex=1:(length(Ptry)); %loop Hyp Exp in sample of 'j', iTry230831n N loop
    parg=Ptry(pindex);Narg=Ntry(Nindex);narg=j;%Arguments for HypExp next
    [HypExpTEMP1(1:(j+1))]=HypExactExp240227Octave(parg,Narg,narg);
    HypExp(Nindex,pindex,(1:j+1))=HypExpTEMP1;%filled to jmax with zeros
    end; % end loop for Hyp Exp
end; % end Loop for Ntry values, calc Hyp Exps for N,n,p

% Bin edges to tally targets-detected, increment 1; -0.5 to +(n+0.5)
    counter1=1;
    for e=-0.5:1:(j+0.5);%"j(Nindex)+1" bins, centres 0...j(Nindex)
    edges(counter1)=e;
    counter1=counter1+1;
    end; % end e loop

%% Jacknife sample
for l=1:L; %sample loci Nit. Inds random in locus*popsize matrix
  %ploc(l)=(sum(indit(l,1:n)))/n;%allele prop locus "l", from sample of n
  for jndx=1:(2*n-j+1);% overlap jac of j,loc'l'
  indalTEMP(l,(1:j))=indalfold(l,(jndx:(jndx+j-1)));
  detNit(l,jndx)=sum(indalTEMP(l,1:j));%TrgtSampd row-loc,col-jack
  end; % end loop for '2*n-j+1' subsamples each sample size j
end; % end loop for sampling over loci 'l'

% Tally results in 'freqdetNit', popsize Nit, & each locus
for l=1:L % loop loci for histogram data.
HistFreq1=histc(detNit(l,:),edges);%in n-j(Nindex)+1 Jacknife reps
HistFreq=HistFreq1(1:(length(HistFreq1)-1));% remove '>max edge' bin
% row vec, bins number detected 0,1,...,(2*n-j+1)
PropNit=HistFreq./sum(HistFreq);
propdetNit(l,(1:j+1))=PropNit;%filled to jmaxALL with zeros
end; % end target-det; dims: locus, bins (0's,1's,...2*n detected)

%% BIAS for each locus, & each hypothesised Ntry, Ptry
for l=1:L; % locus loop, each needs to examine all Ntrys
 for Nindex=1:(length(Ntry)); % loop trial N's each locus,
  for pindex=1:(length(Ptry)); % begin joint estimate p & Nc, loop p values
    % AbsBias, Absolute value:Obs-HypergeoExps, n-1 inds  (count 0,1,...n)
    HypExpTEMP2=(zeros(1,j+1))';
    aa=HypExp(Nindex,pindex,(1:j+1));
    HypExpTEMP2=(squeeze(aa))'; % row vec length j(Nindex)
    PropTEMP=propdetNit(l,:);
    BiasTEMP1=abs(PropTEMP-HypExpTEMP2);%fill to jmax with 0
    BiasTEMP2=sum(BiasTEMP1); % abs bias over all bins, this locus
    BIAS1(Nindex,pindex,l)=BiasTEMP2;% AbsBias all bins, this locus
  end; % end pindex loop
 end; % End loop for simulating and testing Ntry values.
end; %end locus loop.

% Joint Pest-Nest diagnosis each locus, then Nest averaged over all loci
%Diagnose p,N this locus, search min BIAS1, index gives Ptry Ntry indices
for l=1:L;
locBIAS=squeeze(BIAS1(:,:,l)); %bias for this locus (3rd dim)
Min=min(min(locBIAS));
[row,col]=find(locBIAS==Min); % shows row and col of min value
Nest(l)=mean(Ntry(row)); Pest(l)=mean(Ptry(col));% Av if>1 val
end; % end locus loop
% Over loci, AVE and SE of Nest
AveNestDiploid=mean(Nest); SeNestDiploid=(std(Nest))/sqrt(L);
AveNest=AveNestDiploid/2; SeNest=SeNestDiploid/2;% Convert to nbr of inds

%% Now do Nest diagnosis with low CV ""loci"" only, target&(1-target)
SDdetNit=(std(detNit,1,2))';
AvedetNit=(mean(detNit,2))';
CVdetNit=((SDdetNit(:))./(AvedetNit(:)))';
LocNbr(1,1:L)=[1:L];% serialnumbers for ""loci""
CVNest=([CVdetNit(1:L);Nest(1:L);LocNbr(1:L)])';%CVcol1,NestCol2;LocNbrCol3
CVNestSort= (sortrows(CVNest,1))';%sorts by CV, all rows (CV,Nest,LocNbr)
NestLoCV=CVNestSort(2,(1:round(L/10)));%Nest row,lowest 10% of CVs
CVLoLoci=CVNestSort(3,(1:round(L/10)));%""locus"" numbers for AveNestLoCV
% Diploid result, and adjust to numbers of individuals.
AveNestLoCVdiploid=mean(NestLoCV); AveNestLoCV=AveNestLoCVdiploid/2;
SDNestLoCVdiploid=std(NestLoCV); SDNestLoCV=SDNestLoCVdiploid/2;
SENestLoCV=SDNestLoCV/sqrt(L);
DIPminNtry=minNtry/2;DIPmaxNtry=maxNtry/2;% no. of inds, not No. alleles

save ("Output.txt", "DIPminNtry", "DIPmaxNtry", "jj", "AveNestLoCV", "SDNestLoCV", "SENestLoCV");
