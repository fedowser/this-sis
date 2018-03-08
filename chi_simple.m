% n=input('how many data do you want? \n')
n=300;
y=10.^linspace(-2,2,n); % It makes an equal interval in log scale.
%x1 represents chi-prime. And similarly x2 chi-double-prime, xa chi_a, xb
%chi_b.

[x1,x2,xa,xr,x1c,x2c,xac,xrc]=deal(y-y);% I did it to make them dimension-resistant :P. The operator deal is used to assign multiple variable at once.
for ii=1:n
    z=y(ii);
    
    if z<=1
        x1(ii)=z./2-1;
        x2(ii)=2.*z/(3*pi);
        xa(ii)=z./2-1;
        xr(ii)=z./4;
        x1c(ii)=z*(1-5*z/16)-1;
        x2c(ii)=4/3/pi*z*(1-z/2);
        xac(ii)=z-z^2/3-1;
        xrc(ii)=z/2-(z^2)/4;
    else
        x1(ii)=(1/pi).*((z/2-1).*acos(1-2./z)+(-1+4./(3.*z)-4./(3.*z.^2)).*sqrt(z-1));
        x2(ii)=1/3/pi.*(6./z-4./z.^2);
        xa(ii)=-1/2./z;
        
        if z<=2
            xr(ii)=1-z./4-1/2./z;
            xrc(ii)=-1/3/z+1-z/2+z^2/12;
        else
            xr(ii)=1/2./z;
            xrc(ii)=1/3/z;
        end
        
        x1c(ii)=1/pi*((-1+z-5*z^2/16)*acos(1-2/z)+(-19/12+5*z/8+1/z-2/3/z^2)*sqrt(z-1));
        x2c(ii)=4/3/pi/z*(1-1/2/z);
        xac(ii)=-1/3/z;
        
    end
    
end

semilogx(y,x1);
hold on
semilogx(y,x2);
semilogx(y,xa);
semilogx(y,xr);
semilogx(y,x1c,'--');
semilogx(y,x2c,'--');
semilogx(y,xac,'--');
semilogx(y,xrc,'--');
title('AC susceptibility')
hold off
legend('\chi^\prime','\chi^\prime^\prime','\chi_a','\chi_r')
%legend('\chi^\prime, \chi^\prime^\prime','\chi_a, \chi_r','','slab', 'cylinder')
xlabel('$B_a/\mu_0 j_\textnormal{c} R$','Interpreter','latex')
ylabel('$\chi''''$','Interpreter','latex')
%ylabel('\chi^\prime, \chi_a                                                  \chi^\prime^\prime');


%maximum value determination
plotting=[x2;xr;x2c;xrc];
% semilogx(y,plotting)
[ms, ind]=max(plotting,[],2);

% writing on a text file
fileID=fopen('new.txt','w');
fprintf(fileID,'%7s % 9s\r\n','y','maximum');
for ii=ind
fprintf(fileID,'%6.5f   %6.5f\r\n',y(ind),ms);
end
fclose(fileID);
rn={'x2','xr','x2c','xrc'};
cn={'maxval', 'y'}

s=array2table([ms, ind],'RowNames',rn,'VariableNames',cn)
writetable(s,'ne.txt','Delimiter',' ','WriteRowNames',true)
print -dpdf chi.pdf

%https://www.papeeria.com/join?token_id=44bc521a-56ae-4455-8548-d77311064095&retry=3
