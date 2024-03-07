% Hypergeometric Exact Expectations for large N, full freq spectrum 230310
function [HypExp] = HypExactExp230607(parg,Narg,narg); 
% x=numbers detected in narg sampled from Narg with true parg
q=1-parg; % Prop alt allele
qN=round(q*Narg); pN=round(parg*Narg);
% Dim for cumulator
HypExp = zeros(1,(narg+1));  

for x1=1:narg+1; % loop thru bins, +1 so index<>0
  x=x1-1; notx=narg-x; % actual numbers
  notx1=notx+1; %for index, never zero
  %Dimensions for cumulators
  qNpNdetect=ones(1,narg);
  qNundetect=zeros(1,notx);
  pNdetect=zeros(1,x);
  nN=zeros(1,narg);
  nNDenom=ones(1,narg); 

        for det1=1:x; % numerator for Hyp Prob - p part
            det=det1-1; % actual numbers
            pNdetect(det1)=(pN-det)/(x-det);
            if pNdetect(det1)==Inf; pNdetect(det1)=1; end; %Div(0!)trap
        end; % end numerator loop - p part
 
        for undet1=1:notx; % numerator for Hyp Prob - q part
           undet=undet1-1; % actual numbers
             qNundetect(undet1)=(qN-undet)/(notx-undet);
        if qNundetect(undet1)==Inf; qNundetect(undet1)=1; end; %Div(0!)trap
        end; % end numerator loop - q part

  % concat numerator, large to small, traps for all target or non-target
  if notx==narg;
      qNpNdetect(1:notx)=fliplr(qNundetect);
      elseif x==narg;
      qNpNdetect(1:x)=fliplr(pNdetect);
      else
      qNpNdetect(1:notx)=fliplr(qNundetect);
      qNpNdetect((notx+1):narg)=fliplr(pNdetect);
  end; % end concatenation

  for nn=1:narg; %denominator for hyp prob
    nN(nn)=(Narg-nn+1)/(narg-nn+1);
  end; % end denom loop
  nNDenom(1:narg)=fliplr(nN(1:narg)); % large to small, then unity
  
  HypElements=qNpNdetect(:)./nNDenom(:); % column vector, narg rows
  HypExp(x1)=prod(HypElements); %Exp for(Narg,parg,detect)
end; % end main loop
% end function
end