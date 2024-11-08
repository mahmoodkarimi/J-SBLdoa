
function [psaiestofn,nitr]= JSBLprpsd(xmt,scnum,angstpSBL,rfnratio)

% The JSBLprpsd.m contains the Matlab code of the J-SBL doa estimation algoritghm.
%
% psaiestofn: DOA's of sources estimated by the J-SBL algorithm and the fine search
% nitr: Number of iterations in the J-SBL algorithm
% xmt : Matrix containing the outputs of the array
% scnum : Number of sources
% angstpSBL : Anglular step size in the J-SBL algorithm (in degrees)
% rfnratio : Ratio of angstpSBL to the fine search steps


psaimin=-87;% psaimin and psaimax: Minimum and maximum angles of the grid of angles (in degrees)
psaimax=87;
angstpSBLrfn=angstpSBL/rfnratio; % Anglular step size in the grid refinement (in degrees)
[M,ns]=size(xmt); % M: number of sensors, ns: number of snapshots
maxnm=2;%The number of algorithms that exist in this program

psaiestotmp=psaimin+angstpSBL*round((psaimax-psaimin)*(rand(1,scnum))/angstpSBL);%Initial values of psaiesto are randomized 
                                                                                 %integer multiples of angstpSBL
psaiesto=ones(maxnm,1)*psaiestotmp;

nsestd2hato=0;
RestML=(xmt*xmt')/ns; % Estimated covariance matrix
RestMLinv=pinv(RestML);
psaivct=psaimin:angstpSBL:psaimax; % Grid of angles
lnt1=length(psaivct); 
pifcd=1i*pi;
tmp1=pifcd*(sin([psaivct]*pi/180));
pmttst=exp([0:M-1]'*tmp1);% Steering matrix for grid of angles

nitr=500; % Maximum number of iterations
epslnmin=0.001;

%***************************** The SBL algorithm *************************************

gammaold=zeros(lnt1,1);
for jja=1:ns
    gammaold=gammaold+abs(pmttst'*xmt(:,jja));
end
gammaold=gammaold/(M*ns);
Shat=diag(gammaold); %Initial estimate for signal covariance matrix
nsestd2hat=(0.1*(norm(xmt,'fro'))^2)/(M*ns);%Initial estimate for sensors' noise power

epsln=10;
ii=0;
while (epsln>epslnmin)&(ii<nitr)
    ii=ii+1;
    R=pmttst*Shat*pmttst'+nsestd2hat*eye(M);
    Rinv=pinv(R);
    tmp2=RestML*Rinv;
    for jj=1:lnt1
        tmp3=trace(ns*tmp2*pmttst(:,jj)*pmttst(:,jj)'*Rinv)/(ns*pmttst(:,jj)'*Rinv*pmttst(:,jj)+(1/(Shat(jj,jj)+nsestd2hat/M)));
        Shat(jj,jj)=Shat(jj,jj)*tmp3;%Formula of "Multisnapshot Sparse Bayesian ..." paper
                                     % with considering a noninformative prior for Shat(jj,jj)
    end
    Shat=real(Shat);
    gamma=diag(Shat);
    epsln=norm((gamma-gammaold),2)/norm(gammaold,2);
    gammaold=gamma; 
    
    [~,LOCS]= findpeaks(real(diag(Shat)),'SortStr','descend','NPeaks',scnum);
    psaiest=sort(psaivct(LOCS)); % DOA's of sources estimated by SBL method
                                 
    tmp1a=pifcd*(sin([psaiest]*pi/180));
    pmttsta=exp([0:M-1]'*tmp1a);%steering matrix
    tmpr2a=pinv(pmttsta'*pmttsta)*pmttsta';
    nsestd2hat=real(trace((eye(M)-pmttsta*tmpr2a)*RestML)/(M-scnum));%Updating the estimate for noise power
end
figure(40)
spctr=real(diag(Shat));
plot(psaivct,spctr);
grid on
xlabel('Conic angle psai (in degrees)');
str1={'Power of sources estimated by the J-SBL algorithm'};
title(str1);
zoom on
psaiest=sort(psaiest);
psaiestlnt=length(psaiest);
for kk5=1:psaiestlnt
    psaiesto(1,kk5)=psaiest(kk5);
end

nitr=ii;

%************ Applying fine search (Off-grid estimation) to the results of the SBL algorithm **********

R=pmttst*Shat*pmttst'+nsestd2hat*eye(M);
[gammak,LOCS]= findpeaks(real(diag(Shat)),'SortStr','descend','NPeaks',scnum);
psaiest=psaivct(LOCS); % DOA's of sources estimated by SBL method
tmp1a=pifcd*(sin([psaiest]*pi/180));
pmttstb=exp([0:M-1]'*tmp1a);%steering matrix
fngrid=-angstpSBL/2:angstpSBLrfn:angstpSBL/2;
npnts=length(fngrid);

for ja=1:length(psaiest)
    Rmk=R-gammak(ja)*pmttstb(:,ja)*pmttstb(:,ja)';
    Rmkinv=pinv(Rmk);
    tmp2rfn=Rmkinv*RestML*Rmkinv;
    
    psaikvct=psaiest(ja)+fngrid;
    tmp1b=pifcd*(sin([psaikvct]*pi/180));
    pmtfng=exp([0:M-1]'*tmp1b);%steering matrix for the fine grid points 
    for jb=1:npnts
        avctk=pmtfng(:,jb);
        zk=avctk'*Rmkinv*avctk;
        qk=ns*avctk'*tmp2rfn*avctk;
        ggk=2*zk-qk+ns*zk+nsestd2hat*zk*zk*ns/M;
        deltab=ggk*ggk-4*(ns+1)*zk*zk*(1-nsestd2hat*(qk-ns*zk)/M);
        gammakfng(jb)=real((-ggk+sqrt(deltab))/(2*(ns+1)*zk*zk));
        lklhd(jb)=-ns*log(1+zk*gammakfng(jb))+(qk/(zk+(1/gammakfng(jb))))+log(sqrt(ns)/(gammakfng(jb)+nsestd2hat/M));
    end
    [~,ind]=max(real(lklhd));
    psaiestrfn36(ja)=psaikvct(ind);
end
psaiestrfn36=sort(psaiestrfn36);
psaiestrfn36lnt=length(psaiestrfn36);
for kk5=1:psaiestrfn36lnt
    psaiesto(2,kk5)=psaiestrfn36(kk5);
end
psaiestofn=psaiesto(2,:);

end

