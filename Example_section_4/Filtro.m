function [Vq, t0, T0] = Filtro(dp,pc,tc) %

% ENTRADA DE DADOS ----------------------------------------------------------------------------------------------------
%dados=[30 55  tc];%{'Distanciamento entre pontos para Pre-filtragem:','Ponto central para filtragem (Numero impar):',...
% 'Tempo de usinagem estimado no CAM [s]:'};
dados = [dp pc tc]; % [34, 55, 60]
temperature=load('E1.txt');
% SEPARAR O ARQUIVO .txt EM VETORES COLUNA -------------------------------------------------------------------------
termocouple1 = temperature(:,1);
termocouple2 = temperature(:,2);
termocouple3 = temperature(:,3);
termocouple4 = temperature(:,4);
termocouple5 = temperature(:,5);
termocouple6 = temperature(:,6);
termocouple7 = temperature(:,7);
termocouple8 = temperature(:,8);
termocouple9 = temperature(:,9);
termocouple10 = temperature(:,10);
termocouple11 = temperature(:,11);
termocouple12 = temperature(:,12);

%============================================pre-filter: Taking dates with inc=dados(1)
np = length(termocouple1);

inc=dados(1);

cont1=1;
for h=1:inc:np
    termocouple1_pf(cont1)=termocouple1(h);
    termocouple2_pf(cont1)=termocouple2(h);
    termocouple3_pf(cont1)=termocouple3(h);
    termocouple4_pf(cont1)=termocouple4(h);
    termocouple5_pf(cont1)=termocouple5(h);
    termocouple6_pf(cont1)=termocouple6(h);
    termocouple7_pf(cont1)=termocouple7(h);
    termocouple8_pf(cont1)=termocouple8(h);
    termocouple9_pf(cont1)=termocouple9(h);
    termocouple10_pf(cont1)=termocouple10(h);
    termocouple11_pf(cont1)=termocouple11(h);
    termocouple12_pf(cont1)=termocouple12(h);
    cont1=cont1+1;
end


%=======================================Filter===========================================================================
nppf=length(termocouple1_pf);

vu=dados(2)-1;
vt=2*vu+1;


for i=dados(2):nppf-vu
    soma1=0;
    soma2=0;
    soma3=0;
   soma4=0;
    soma5=0;
    soma6=0;
     soma7=0;
    soma8=0;
    soma9=0;
     soma10=0;
    soma11=0;
    soma12=0;
    for j=i-vu:i+vu
        soma1=soma1+termocouple1_pf(j);
        soma2=soma2+termocouple2_pf(j);
        soma3=soma3+termocouple3_pf(j);
         soma4=soma4+termocouple4_pf(j);
        soma5=soma5+termocouple5_pf(j);
        soma6=soma6+termocouple6_pf(j);
         soma7=soma7+termocouple7_pf(j);
        soma8=soma8+termocouple8_pf(j);
        soma9=soma9+termocouple9_pf(j);
         soma10=soma10+termocouple10_pf(j);
        soma11=soma11+termocouple11_pf(j);
        soma12=soma12+termocouple12_pf(j);
        
    end
    termocouple1_f(i)=soma1/vt;
    termocouple2_f(i)=soma2/vt;
    termocouple3_f(i)=soma3/vt;
    termocouple4_f(i)=soma4/vt;
    termocouple5_f(i)=soma5/vt;
    termocouple6_f(i)=soma6/vt;
   termocouple7_f(i)=soma7/vt;
    termocouple8_f(i)=soma8/vt;
    termocouple9_f(i)=soma9/vt;
   termocouple10_f(i)=soma10/vt;
    termocouple11_f(i)=soma11/vt;
    termocouple12_f(i)=soma12/vt;
    
end


% [vmax1,imax1]=max(termocouple1_f);
% [vmax2,imax2]=max(termocouple2_f);
% [vmax3,imax3]=max(termocouple3_f);
% [vmax4,imax4]=max(termocouple4_f);
% [vmax5,imax5]=max(termocouple5_f);
% [vmax6,imax6]=max(termocouple6_f);
% [vmax7,imax7]=max(termocouple7_f);
% [vmax8,imax8]=max(termocouple8_f);
% [vmax9,imax9]=max(termocouple9_f);
% [vmax10,imax10]=max(termocouple10_f);
% [vmax11,imax11]=max(termocouple11_f);
[vmax12,imax12]=max(termocouple12_f);
[vmin1,imin1]=min(termocouple1_f(dados(2):length(termocouple1_f)));
%========== Time for max temperature in termocouple 12

npm=length((termocouple1_f(dados(2):imax12)));
         
inc_time_m=dados(3)/(npm-1);

time_m(1)=0;


for y=1:npm-1
    time_m(y+1)=y*inc_time_m;
    
end
t0=time_m;
T0=vmin1;
Vq=[termocouple1_f(dados(2):imax12); termocouple2_f(dados(2):imax12); termocouple3_f(dados(2):imax12); termocouple4_f(dados(2):imax12); termocouple5_f(dados(2):imax12); termocouple6_f(dados(2):imax12); termocouple7_f(dados(2):imax12); termocouple8_f(dados(2):imax12); termocouple9_f(dados(2):imax12); termocouple10_f(dados(2):imax12); termocouple11_f(dados(2):imax12); termocouple12_f(dados(2):imax12)];

end

% figure('NumberTitle','off','Name','Filtered Temperature curves','Position',[80 80 850 530]);
% plot(time_m,termocouple1_f(dados(2):imax12));
% hold on
% plot(time_m,termocouple2_f(dados(2):imax12));
% hold on
% plot(time_m,termocouple3_f(dados(2):imax12));
% hold on
% plot(time_m,termocouple4_f(dados(2):imax12));
% hold on
% plot(time_m,termocouple5_f(dados(2):imax12));
% hold on
% plot(time_m,termocouple6_f(dados(2):imax12));
% hold on
% plot(time_m,termocouple7_f(dados(2):imax12));
% hold on
% plot(time_m,termocouple8_f(dados(2):imax12));
% hold on
% plot(time_m,termocouple9_f(dados(2):imax12));
% hold on
% plot(time_m,termocouple10_f(dados(2):imax12));
% hold on
% plot(time_m,termocouple11_f(dados(2):imax12));
% hold on
% plot(time_m,termocouple12_f(dados(2):imax12));
% title('Curvas filtradas de aquecimento dos corpos-de-prova ate a amplitude maxima no termopar 12')
% xlabel('Time')
% ylabel('Temperatura [C]')
% legend('termopar 1','termopar 2','termopar 3','termopar 4','termopar 5','termopar 6', 'termopar 7','termopar 9','termopar 10', 'termopar 11','termopar 12')
% zoom on


 
%==============================================================



