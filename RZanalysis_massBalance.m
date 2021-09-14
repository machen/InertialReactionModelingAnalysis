clear all; 
clc;
close all

caseDir = [pwd '/'];

saveDir = [caseDir 'postProcess/'];


%% read constant/transportProperties, specifically the viscosity
fidNu = fopen([caseDir 'constant/transportProperties'],'r'); 
g = textscan(fidNu, '%s','delimiter','\n'); fclose(fidNu);
g = g{1};

iNu = []; iValue = [];
for iL = 1:length(g)
    if contains(g{iL}, 'nu')
        iNu= iL;
    end
end

tmp = textscan(char(g{iNu} ), '%s','delimiter',']'); tmp = tmp{1}(2);
tmp = textscan(char(tmp ), '%s','delimiter',';'); 
nu_phy = str2num(char(tmp{1})); % [m^2/s] 
                                % 0.8926*1e-6: physical water kinematic ...
                                %           viscosity at 25 celsius degree

clear g

%% ReIn, ID the Re number
kRe = strfind(caseDir,'Re');
ReIn = str2double(caseDir(kRe+2:end-1));

%% load cell velocity, loads the velocity components (UxC and UyC)
load([saveDir 'UC.mat'],'UxC','UyC');

Unorm = sqrt(UxC.*UxC + UyC.*UyC); % Magnitude of velocity.

%% H & geometry
load([saveDir 'coordinates.mat']);%,'ndx','ny','nz','Lx','Ly','Lz',...
%                                'dx','dy','dz','xFaces','yFaces','zFaces',...
%                                'xCells','yCells','zCells',...
%                                'nPoints','nCells','nFaces','nInternalFaces',...
%                                'convertToMeters');

b0 = sum(~isnan(UxC(:,1)))*dx; %[m] aperture || mean separation    
                   
kH = strfind(caseDir, 'H');
knx = strfind(caseDir,'_nx');
H = str2double(caseDir(kH+1:knx-1));
    

    
%% RZ analysis

    fidQio = fopen([caseDir '/postProcess/Qin_outEvolv.dat'], 'w'); 
       
    lamda=zeros(length(xCells),2); %cell(iH,iRe)
    xBaseAll=zeros(length(xCells),2); yBaseAll=zeros(length(xCells),2);
    xBaseMain=zeros(length(xCells),2); yBaseMain=zeros(length(xCells),2);
    
    ixRZ =zeros(sum(~isnan(UxC(:))), 1);
    iyRZ =ixRZ;
    lengthRZ = 0;
    
    for ix = 1:length(xCells)        
        %% ============== apparent aperture measure
        indPore = find(~isnan(UxC(:,ix)));

        Uflow = [0; UxC(indPore,ix); 0];

        xBase = ones(length(indPore)+2,1)*xCells(ix);
        yBase = [yCells(indPore(1))-dy/2; ...
                 yCells(indPore); 
                 yCells(indPore(end))+dy/2];

        baseCoord = yBase-yBase(1);

        lamda(ix,1) = baseCoord(end);            

        baseCoordIntrp = linspace(baseCoord(1),baseCoord(end),10000);
        UflowIntrp = interpn(baseCoord,Uflow,baseCoordIntrp','linear');            

        cumQempirical1 = cumtrapz(baseCoordIntrp,UflowIntrp);
        indMIMintfce1 = find([0; abs(diff(sign(cumQempirical1)))]==2);

        cumQempirical2 = cumtrapz(baseCoordIntrp,flip(UflowIntrp));
        indMIMintfce2 = find([0; abs(diff(sign(cumQempirical2)))]==2);

        if isempty(indMIMintfce1);  indMIMintfce1 =0; end
        if isempty(indMIMintfce2);  indMIMintfce2 =0; end

        indMIMintfce1 = indMIMintfce1(end);
        indMIMintfce2 = indMIMintfce2(end);

        Qempirical = trapz(baseCoordIntrp(indMIMintfce1+1:end-indMIMintfce2),...
                               UflowIntrp(indMIMintfce1+1:end-indMIMintfce2));            
        
        if ix == 1
            Qinlet = Qempirical; 
        elseif ix == length(xCells)
            Qoutlet = Qempirical; 
        end
        Re = Qempirical/nu_phy;
        
                           
        lamda(ix,2) = lamda(ix,1)*...
            (length(UflowIntrp)-indMIMintfce1-indMIMintfce2)/length(UflowIntrp);

        xBaseMain(ix+1,:) = ...
           [xBase(1)+(xBase(end)-xBase(1))*indMIMintfce1/length(UflowIntrp),...
            xBase(end)-(xBase(end)-xBase(1))*indMIMintfce2/length(UflowIntrp)];
        yBaseMain(ix+1,:) = ...
            [yBase(1)+(yBase(end)-yBase(1))*indMIMintfce1/length(UflowIntrp),...
             yBase(end)-(yBase(end)-yBase(1))*indMIMintfce2/length(UflowIntrp)];

        xBaseAll(ix+1,:) = [xBase(1), xBase(end)];
        yBaseAll(ix+1,:) = [yBase(1), yBase(end)];
%         
%             clf;
%             xverts = [baseCoordIntrp(1:end-1); baseCoordIntrp(1:end-1); baseCoordIntrp(2:end); baseCoordIntrp(2:end)];
%             yverts = [zeros(1,numel(UflowIntrp)-1); UflowIntrp(1:end-1)'; UflowIntrp(2:end)'; zeros(1,numel(UflowIntrp)-1)];
%             p = patch(xverts,yverts,'b','LineWidth',1.5);
%             xlim([0,0.001])
%             drawnow; pause(0.01);

%         xCenter = 2\(xBase(1)+xBase(end));
%         yCenter = 2\(yBase(1)+yBase(end));
%         sqD = (xBase-xCenter).^2 + (yBase-yCenter).^2;
% 
%         flowQTheory = Qinlet;%UmeanPoisTheory * b0;            
%         UflowTheory = 6*flowQTheory/ lamda(ix,1)^3 * ...
%             ((lamda(ix,1)/2)^2 - sqD);
%         UflowIntrpTheory = interpn(baseCoord,UflowTheory,baseCoordIntrp','linear');
%         Qpredict = trapz(baseCoord,UflowTheory);
% 
%         Uflow = Uflow/Umean*dx*figSize/10;  % rescaled profile            
%         xPara = xBase + Uflow;
%         yPara = yBase;
% 
%         UflowTheory = UflowTheory/Umean*dx*figSize/10;% rescaled profile             
%         xParaTheory = xBase + UflowTheory;
%         yParaTheory = yBase;


        txt = sprintf('x:%2.5f | Qout/in:%.4f | H:%.1f ReIn:%3.3f | Re:%2.4f \n',...
                xCells(ix), Qempirical/Qinlet, H, ReIn, Re);
        fprintf(fidQio,txt);

        clear Uflow;

        %% ============== Figure
        
        if sum(diff([xBaseAll(ix,1); xBaseMain(ix,1)])) ||sum(diff([yBaseAll(ix,1); yBaseMain(ix,1)]))
            
            indRZ = indPore(yCells(indPore)<=yBaseMain(ix,1));
            if ~isempty(indRZ)
                iyRZ(lengthRZ+1:lengthRZ+length(indRZ)) = indRZ;
                ixRZ(lengthRZ+1:lengthRZ+length(indRZ)) = ix;
                lengthRZ = lengthRZ + length(indRZ);
            end    
        end

        if sum(diff([xBaseAll(ix,2); xBaseMain(ix,2)])) ||sum(diff([yBaseAll(ix,2); yBaseMain(ix,2)]))
            indRZ = indPore(yCells(indPore)>=yBaseMain(ix,2));
            if ~isempty(indRZ)
                iyRZ(lengthRZ+1:lengthRZ+length(indRZ)) = indRZ;
                ixRZ(lengthRZ+1:lengthRZ+length(indRZ)) = ix;
                lengthRZ = lengthRZ + length(indRZ);
            end    
        end
        
    end
    
    Areas = sum(lamda);
    RZratio =  abs(diff(Areas))/max(Areas);

    ixRZ = ixRZ(ixRZ~=0);
    iyRZ = iyRZ(iyRZ~=0);
    
    Qin_out = Qoutlet / Qinlet;
    
    RZind = (ixRZ-1)*ny+iyRZ;
    UmeanRZ = mean(Unorm(RZind));
    b0RZ = b0 * sqrt(RZratio);
    
    ReRZ = UmeanRZ*b0RZ/nu_phy;

    
    
    ttxt  = sprintf('H: %.1f  Re: %.2f  (out/in: %.4f  Re: %.5f   RZ/Pore: %.5f) \n',...
        H, ReIn, Qin_out, Re, RZratio);
    
    fprintf(fidQio,[ttxt]); 
    fclose(fidQio);
  
    save([caseDir '/postProcess/RZresults_massBalance.mat'],'lamda','xBaseAll','yBaseAll',...
                                                  'xBaseMain','yBaseMain',...
                                                  'ixRZ','iyRZ','ReRZ',...
                                                  'Qinlet',...
                                                  'Re','Qin_out','RZratio');
    
                                              




%% Quit MATLAB
quit;