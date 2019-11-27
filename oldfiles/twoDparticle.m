clear;
clc;
tic;
for ikpj=10001
    close all;
    %% Initial Parameter 
    idum0=ikpj;
    fileName=['Part2D_' num2str(ikpj)] ;
    addInterface=false;
    YLmemb=true;
    
    pVolume=0.4;                    % volume content 0.25
    nPart=100;
    nTry=10000;
    
    xBound=10.0;
    yBound=xBound;
%     zBound=1.0;
    
    TOLERDIST=0.1;
    RDISPERSE=0;                    % dimension disperse (if radium have no difference, set 0)
    
    NMAX=1000;
    Nseed0=24;
    Nseed=24;
    
    idum1=idum0+1;                  %
    idum2=idum1+1;                  % The random
    idum3=idum2+1;                  % seeds
    idum4=idum3+1;                  %
    idum5=idum4+1;                  %
    
    radius=sqrt(4.0*xBound*yBound*pVolume/nPart/pi);
    RINTER=2*pi*radius/(Nseed0);
    
    dirname=['nPart_',num2str(nPart),'_pVolume_',num2str(pVolume)];
    if ~exist(dirname)
        mkdir(dirname);
    end
    
    lengElem0=2.0*pi*radius/Nseed;
    lengElem1=lengElem0*2.0;
    lengElem2=lengElem0*sqrt(2);
    
    
    % modify radius
    myStream=RandStream('mt19937ar','Seed',idum0);
    r0=rand(myStream,1,nPart);
    Rdi(1:nPart)=radius+RDISPERSE*abs(r0(1:nPart));
    Rdi2(1:nPart)=Rdi(1:nPart).^2;
    TEMP=sum(Rdi2);
    TEMP=sqrt(TEMP/nPart/radius^2);
    Rdi(1:nPart)=Rdi(1:nPart)/TEMP;
    rmin=min(Rdi)*RINTER/10;
%     xyAlready(1,1:5)=[0 0 2.0^2*xBound 2.0^2*yBound radius];
    
    Nleng=nPart*nTry;
    myStream=RandStream('mt19937ar','Seed',idum1);
    xRand=rand(myStream,1,Nleng);
    myStream=RandStream('mt19937ar','Seed',idum2);
    yRand=rand(myStream,1,Nleng);
    imn=1;
    thisGen=true;
    
    
    % initial arrays and cell
    xyAlreadyCell=cell(nPart,1);
    xCrossArray=[];
    yCrossArray=[];
    cornerCrossArray=[];
    
    figure(1);
    %% Fiber Part
    for kpart=1:nPart
%         part=kpart;
        if thisGen==false
            warning('this seed is not good for this simulation');
            break
        end
        for knum=1:nTry

            nGen=false;
            
            % random generation
            xT=xBound*(2.0*xRand(imn)-1.0);
            yT=yBound*(2.0*yRand(imn)-1.0);
            
            imn=imn+1;
%             deltaR=2.0*Rdi(kpart)*TOLERDIST;
            deltaR=2.0*Rdi(kpart)*TOLERDIST;
            r1=Rdi(kpart)-deltaR;
            r2=Rdi(kpart)*(1+RINTER)+deltaR;
            
            if xT>xBound || xT<-xBound  % x+
                continue;
            end
            if yT>yBound || yT<-yBound  % x+
                continue;
            end
            
            if xBound-xT>r1 && xBound-xT<r2  % x+
                continue;
            end
            
            if yBound-yT>r1 && yBound-yT<r2  % y+
                continue;
            end
            
            if xBound+xT>r1 && xBound+xT<r2  % x-
                continue;
            end
            
            if yBound+yT>r1 && yBound+yT<r2  % y-
                continue;
            end
            
            if xBound-xT<Rdi(kpart)
                if yBound-yT<Rdi(kpart)
                    xyTemp(1,1:4)=[1 xT-2.0*xBound yT-2.0*yBound Rdi(kpart)];
                    xyTemp(2,1:4)=[2 xT yT-2.0*yBound Rdi(kpart)];
                    xyTemp(3,1:4)=[3 xT yT Rdi(kpart) ];
                    xyTemp(4,1:4)=[4 xT-2.0*xBound yT Rdi(kpart)];
                    crossType='Corner';
                elseif yBound+yT<Rdi(kpart)
                    xyTemp(1,1:4)=[1 xT-2.0*xBound yT Rdi(kpart)];
                    xyTemp(2,1:4)=[2 xT yT Rdi(kpart)];
                    xyTemp(3,1:4)=[3 xT yT+2.0*yBound Rdi(kpart)];
                    xyTemp(4,1:4)=[4 xT-2.0*xBound yT+2.0*yBound Rdi(kpart)];
                    crossType='Corner';
                else
                    xyTemp(1,1:4)=[1 xT-2.0*xBound yT Rdi(kpart)];
                    xyTemp(2,1:4)=[2 xT yT Rdi(kpart)];
                    crossType='X';
                end
                
            elseif xBound+xT<Rdi(kpart)
                if yBound-yT<Rdi(kpart)
                    xyTemp(1,1:4)=[1 xT yT-2.0*yBound Rdi(kpart)];
                    xyTemp(2,1:4)=[2 xT+2.0*xBound yT-2.0*yBound Rdi(kpart)];
                    xyTemp(3,1:4)=[3 xT+2.0*xBound yT Rdi(kpart)];
                    xyTemp(4,1:4)=[4 xT yT Rdi(kpart)];
                    crossType='Corner';
                elseif yBound+yT<Rdi(kpart)
                    xyTemp(1,1:4)=[1 xT yT Rdi(kpart)];
                    xyTemp(2,1:4)=[2 xT+2.0*xBound yT Rdi(kpart)];
                    xyTemp(3,1:4)=[3 xT+2.0*xBound yT+2.0*yBound Rdi(kpart)];
                    xyTemp(4,1:4)=[4 xT yT+2.0*yBound Rdi(kpart)];
                    crossType='Corner';
                else
                    xyTemp(1,1:4)=[1 xT yT Rdi(kpart)];
                    xyTemp(2,1:4)=[2 xT+2.0*xBound yT Rdi(kpart)];
                    crossType='X';
                end
            else
                if yBound-yT<Rdi(kpart)
                    xyTemp(1,1:4)=[1 xT yT-2.0*yBound Rdi(kpart)];
                    xyTemp(2,1:4)=[2 xT yT Rdi(kpart)];
                    crossType='Y';
                elseif yBound+yT<Rdi(kpart)
                    xyTemp(1,1:4)=[1 xT yT Rdi(kpart)];
                    xyTemp(2,1:4)=[2 xT yT+2.0*yBound Rdi(kpart)];
                    crossType='Y';
                else
                    xyTemp(1,1:4)=[1 xT yT Rdi(kpart)];
                    crossType='None';
                    %kpart
                end
            end
            
            xyAlready=[];
            for lpart=1:kpart-1
                xyAlready=[xyAlready;xyAlreadyCell{lpart}];
            end
            leng1=size(xyTemp,1);
            leng2=size(xyAlready,1);
            for nTemp=1:leng1                                           % check if contact other fibers
                x1=xyTemp(nTemp,2);
                y1=xyTemp(nTemp,3);
                r1=xyTemp(nTemp,4);
                for kTemp=1:leng2
                    x2=xyAlready(kTemp,1);
                    y2=xyAlready(kTemp,2);
                    r2=xyAlready(kTemp,3);
                    dist=sqrt((x1-x2)^2+(y1-y2)^2);
                    if dist<(r1+r2)*(1.0+TOLERDIST)
                        nGen=true;
                        break
                    end
                end
                if nGen==true
                    break
                end
            end
            if nGen==true  % x+   
                clear xyTemp;
                continue;
            else
%                 xyAlready(leng2+1:leng2+leng1,3:5)=xyTemp(1:leng1,2:4);
%                 xyAlready(leng2+1:leng2+leng1,2)=kpart;
                xyAlreadyCell{kpart}=xyTemp(1:leng1,2:4);
                for kp=1:leng1
                    text(xyTemp(kp,2),xyTemp(kp,3),num2str(kpart)); hold on
                    plot_circle(xyTemp(kp,2),xyTemp(kp,3), Rdi(kpart), 'r'); hold on
                end
                
                switch crossType
                    case 'Corner'
                        cornerCrossArray=[cornerCrossArray,kpart];
                    case 'X'
                        xCrossArray=[xCrossArray,kpart];
                    case 'Y'
                        yCrossArray=[yCrossArray,kpart];
                    case 'None'  
                end
                
                clear xyTemp xyAlready;
                break;
            end
   
        end
        if (knum==nTry)&&(nGen==true)
%             pause
            thisGen=false;
        end
    end
    if thisGen==false
        clear all;
        continue;
    end
    xtemp=[-xBound xBound xBound -xBound -xBound];
    ytemp=[-yBound -yBound yBound yBound -yBound];
    plot(xtemp,ytemp,'g');  hold on;
    
    EdgeX=[];
    EdgeY=[];
    CellEdgePoint=cell(nPart,1);
    % Edge X:
    if length(cornerCrossArray)==1
        x0temp=xyAlreadyCell{cornerCrossArray}(4,1);
        y0temp=xyAlreadyCell{cornerCrossArray}(4,2);
        r0temp=xyAlreadyCell{cornerCrossArray}(4,3);
        EdgeX(1)=y0temp-sqrt(r0temp^2-(-xBound-x0temp)^2);
        CellEdgePoint{cornerCrossArray}=[CellEdgePoint{cornerCrossArray};[-xBound,EdgeX(1)]];
        clear x0temp y0temp r0temp
    end
    for kedge=1:length(xCrossArray)
        x0temp=xyAlreadyCell{xCrossArray(kedge)}(1,1);
        y0temp=xyAlreadyCell{xCrossArray(kedge)}(1,2);
        r0temp=xyAlreadyCell{xCrossArray(kedge)}(1,3);
        EdgeX=[EdgeX,y0temp+sqrt(r0temp^2-(-xBound-x0temp)^2),y0temp-sqrt(r0temp^2-(-xBound-x0temp)^2)];
        CellEdgePoint{xCrossArray(kedge)}=[CellEdgePoint{xCrossArray(kedge)};[-xBound,EdgeX(end-1)];[-xBound,EdgeX(end)];];
        clear x0temp y0temp r0temp
    end
    if length(cornerCrossArray)==1
        x0temp=xyAlreadyCell{cornerCrossArray}(1,1);
        y0temp=xyAlreadyCell{cornerCrossArray}(1,2);
        r0temp=xyAlreadyCell{cornerCrossArray}(1,3);
        EdgeX=[EdgeX,y0temp+sqrt(r0temp^2-(-xBound-x0temp)^2)];
        CellEdgePoint{cornerCrossArray}=[CellEdgePoint{cornerCrossArray};[-xBound,EdgeX(end)]];
        clear x0temp y0temp r0temp
    end
    
    % Edge Y:
    if length(cornerCrossArray)==1
        x0temp=xyAlreadyCell{cornerCrossArray}(4,1);
        y0temp=xyAlreadyCell{cornerCrossArray}(4,2);
        r0temp=xyAlreadyCell{cornerCrossArray}(4,3);
        EdgeY(1)=x0temp+sqrt(r0temp^2-(yBound-y0temp)^2);
        CellEdgePoint{cornerCrossArray}=[CellEdgePoint{cornerCrossArray};[EdgeY(1),yBound]];
        clear x0temp y0temp r0temp
    end
    for kedge=1:length(yCrossArray)
        x0temp=xyAlreadyCell{yCrossArray(kedge)}(2,1);
        y0temp=xyAlreadyCell{yCrossArray(kedge)}(2,2);
        r0temp=xyAlreadyCell{yCrossArray(kedge)}(2,3);
        EdgeY=[EdgeY,x0temp-sqrt(r0temp^2-(yBound-y0temp)^2),x0temp+sqrt(r0temp^2-(yBound-y0temp)^2)];
        CellEdgePoint{yCrossArray(kedge)}=[CellEdgePoint{yCrossArray(kedge)};[EdgeY(end-1),yBound];[EdgeY(end),yBound];];
        clear x0temp y0temp r0temp
    end
    if length(cornerCrossArray)==1
        x0temp=xyAlreadyCell{cornerCrossArray}(3,1);
        y0temp=xyAlreadyCell{cornerCrossArray}(3,2);
        r0temp=xyAlreadyCell{cornerCrossArray}(3,3);
        EdgeY=[EdgeY,x0temp-sqrt(r0temp^2-(yBound-y0temp)^2)];
        CellEdgePoint{cornerCrossArray}=[CellEdgePoint{cornerCrossArray};[EdgeY(end),yBound]];
        clear x0temp y0temp r0temp
    end
    EdgeX=sort(EdgeX,'descend');
    plot(-xBound*ones(1,length(EdgeX)),EdgeX,'LineStyle','none','Marker','o','MarkerEdgeColor','b'); hold on;
    plot(xBound*ones(1,length(EdgeX)),EdgeX,'LineStyle','none','Marker','o','MarkerEdgeColor','b'); hold on;
    plot(EdgeY,yBound*ones(1,length(EdgeY)),'LineStyle','none','Marker','o','MarkerEdgeColor','k'); hold on;
    plot(EdgeY,-yBound*ones(1,length(EdgeY)),'LineStyle','none','Marker','o','MarkerEdgeColor','k'); hold on;
    EdgeY=sort(EdgeY);
    EdgeXmiddle=[xBound,EdgeX]*0.5+[EdgeX,-xBound]*0.5;
    EdgeYmiddle=[-yBound,EdgeY]*0.5+[EdgeY,yBound]*0.5;
    plot(-xBound*ones(1,length(EdgeXmiddle)),EdgeXmiddle,'LineStyle','none','Marker','x','MarkerEdgeColor','b'); hold on;
    plot(xBound*ones(1,length(EdgeXmiddle)),EdgeXmiddle,'LineStyle','none','Marker','x','MarkerEdgeColor','b'); hold on;
    plot(EdgeYmiddle,yBound*ones(1,length(EdgeYmiddle)),'LineStyle','none','Marker','x','MarkerEdgeColor','k'); hold on;
    plot(EdgeYmiddle,-yBound*ones(1,length(EdgeYmiddle)),'LineStyle','none','Marker','x','MarkerEdgeColor','k'); hold on;
%     hold off;
%     grid on;
%     axis equal 
    
    %% CNF part
    
    %% Python Script Writer
    tolAba=1e-6; % minimine float in Abaqus/CAE
    fID=fopen([dirname,'\',fileName,'.py'],'w');
    fprintf(fID,'%s \n',['# -*- coding: mbcs -*-']);
    fprintf(fID,'%s \n',['from abaqus import *']);
    fprintf(fID,'%s \n',['from abaqusConstants import *']);
    fprintf(fID,'%s \n',['from caeModules import *']);
    fprintf(fID,'%s \n',['from driverUtils import executeOnCaeStartup']);
    fprintf(fID,'%s \n',['from textRepr import *']);
    fprintf(fID,'%s \n',['# Initial Parameter']);
    fprintf(fID,'%s \n',['xBound=',num2str(xBound)]);
    fprintf(fID,'%s \n',['yBound=',num2str(yBound)]);
    fprintf(fID,'%s \n',['Mdb()']);
    fprintf(fID,'%s \n',['']);
    % Create 2D part
    fprintf(fID,'%s \n',['# Create 2D Part']);
    fprintf(fID,'%s \n',['s = mdb.models[''Model-1''].ConstrainedSketch(name=''__profile__'',sheetSize=',num2str(20*xBound),')']);
    fprintf(fID,'%s \n',['s.setPrimaryObject(option=STANDALONE)']);
    fprintf(fID,'%s \n',['s.rectangle(point1=(-xBound, -yBound), point2=(xBound,yBound))']);
    fprintf(fID,'%s \n',['p = mdb.models[''Model-1''].Part(name=''Part-1'', dimensionality=TWO_D_PLANAR,type=DEFORMABLE_BODY)']);
    fprintf(fID,'%s \n',['p = mdb.models[''Model-1''].parts[''Part-1'']']);
    fprintf(fID,'%s \n',['p.BaseShell(sketch=s)']);
    fprintf(fID,'%s \n',['s.unsetPrimaryObject()']);
    fprintf(fID,'%s \n',['p = mdb.models[''Model-1''].parts[''Part-1'']']);
    fprintf(fID,'%s \n',['session.viewports[''Viewport: 1''].setValues(displayedObject=p)']);
    fprintf(fID,'%s \n',['del mdb.models[''Model-1''].sketches[''__profile__'']']);
    fprintf(fID,'%s \n',['']);
    fprintf(fID,'%s \n',['# Partition Plane']);
    fprintf(fID,'%s \n',['p = mdb.models[''Model-1''].parts[''Part-1'']']);
    fprintf(fID,'%s \n',['f, e, d1 = p.faces, p.edges, p.datums']);
    fprintf(fID,'%s \n',['t = p.MakeSketchTransform(sketchPlane=f[0], sketchPlaneSide=SIDE1, origin=(0.0, 0.0, 0.0))']);
    fprintf(fID,'%s \n',['s1 = mdb.models[''Model-1''].ConstrainedSketch(name=''__profile__'', sheetSize=',num2str(20*xBound),', gridSpacing=',num2str(xBound),', transform=t)']);
    fprintf(fID,'%s \n',['g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints']);
    fprintf(fID,'%s \n',['s1.setPrimaryObject(option=SUPERIMPOSE)']);
    fprintf(fID,'%s \n',['p = mdb.models[''Model-1''].parts[''Part-1'']']);
    fprintf(fID,'%s \n',['p.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)']);
    % Draw the Particles
    fprintf(fID,'%s \n',['# Draw the Particles']);
    InterfacePointMatrix=[];  % Matrix recording points on the interfaces(edges)
    for kpart=1:nPart
        fprintf(fID,'%s \n',['# Part ',num2str(kpart)]);
        if ~ismember(kpart,[cornerCrossArray,xCrossArray,yCrossArray])  % inner fiber
            x0temp=num2str(xyAlreadyCell{kpart}(1,1));
            y0temp=num2str(xyAlreadyCell{kpart}(1,2));
            x1temp=num2str(xyAlreadyCell{kpart}(1,1));
            y1temp=num2str(xyAlreadyCell{kpart}(1,2)-xyAlreadyCell{kpart}(1,3));
            InterfacePointMatrix=[InterfacePointMatrix;[xyAlreadyCell{kpart}(1,1),xyAlreadyCell{kpart}(1,2)+xyAlreadyCell{kpart}(1,3)]];
            fprintf(fID,'%s \n',['s1.ArcByCenterEnds(center=(',x0temp,', ',y0temp,'), point1=(',x1temp,', ',y1temp,'), point2=(',x1temp,', ',y1temp,'), direction=COUNTERCLOCKWISE)']);
        elseif ismember(kpart,cornerCrossArray) % Corner fiber   
            x0temp=num2str(xyAlreadyCell{kpart}(1,1));
            y0temp=num2str(xyAlreadyCell{kpart}(1,2));
            r0temp=xyAlreadyCell{kpart}(1,3);
            vectortemp=[-xBound;-yBound]-[xyAlreadyCell{kpart}(1,1);xyAlreadyCell{kpart}(1,2)];
            vectortemp=(dot(vectortemp,[-1;-1])/abs(dot(vectortemp,[-1;-1]))*vectortemp/norm(vectortemp));
            x1temp=num2str(xyAlreadyCell{kpart}(1,1)+r0temp*vectortemp(1));
            y1temp=num2str(xyAlreadyCell{kpart}(1,2)+r0temp*vectortemp(2));
            fprintf(fID,'%s \n',['s1.ArcByCenterEnds(center=(',x0temp,', ',y0temp,'), point1=(',x1temp,', ',y1temp,'), point2=(',x1temp,', ',y1temp,'), direction=COUNTERCLOCKWISE)']);
            clear x0temp y0temp x1temp y1temp vectortemp 
            x0temp=xyAlreadyCell{kpart}(1,1);
            y0temp=xyAlreadyCell{kpart}(1,2);
            x2temp=CellEdgePoint{kpart}(2,1);
            y2temp=CellEdgePoint{kpart}(2,2);
            x3temp=CellEdgePoint{kpart}(3,1);
            y3temp=CellEdgePoint{kpart}(3,2)-2*yBound;
            vectortemp=0.5*[x2temp-x0temp;y2temp-y0temp]+0.5*[x3temp-x0temp;y3temp-y0temp];
            vectortemp=vectortemp/norm(vectortemp);
            InterfacePointMatrix=[InterfacePointMatrix;[x0temp+vectortemp(1)*r0temp,y0temp+vectortemp(2)*r0temp]];
            clear x0temp y0temp x2temp y2temp x3temp y3temp vectortemp r0temp
            
            x0temp=num2str(xyAlreadyCell{kpart}(2,1));
            y0temp=num2str(xyAlreadyCell{kpart}(2,2));
            r0temp=xyAlreadyCell{kpart}(2,3);
            vectortemp=[xBound;-yBound]-[xyAlreadyCell{kpart}(2,1);xyAlreadyCell{kpart}(2,2)];
            vectortemp=(dot(vectortemp,[1;-1])/abs(dot(vectortemp,[1;-1]))*vectortemp/norm(vectortemp));
            x1temp=num2str(xyAlreadyCell{kpart}(2,1)+r0temp*vectortemp(1));
            y1temp=num2str(xyAlreadyCell{kpart}(2,2)+r0temp*vectortemp(2));
            fprintf(fID,'%s \n',['s1.ArcByCenterEnds(center=(',x0temp,', ',y0temp,'), point1=(',x1temp,', ',y1temp,'), point2=(',x1temp,', ',y1temp,'), direction=COUNTERCLOCKWISE)']);
            clear x0temp y0temp x1temp y1temp vectortemp 
            x0temp=xyAlreadyCell{kpart}(2,1);
            y0temp=xyAlreadyCell{kpart}(2,2);
            x2temp=CellEdgePoint{kpart}(4,1);
            y2temp=CellEdgePoint{kpart}(4,2)-2*yBound;
            x3temp=CellEdgePoint{kpart}(2,1)+2*xBound;
            y3temp=CellEdgePoint{kpart}(2,2);
            vectortemp=0.5*[x2temp-x0temp;y2temp-y0temp]+0.5*[x3temp-x0temp;y3temp-y0temp];
            vectortemp=vectortemp/norm(vectortemp);
            InterfacePointMatrix=[InterfacePointMatrix;[x0temp+vectortemp(1)*r0temp,y0temp+vectortemp(2)*r0temp]];
            clear x0temp y0temp x2temp y2temp x3temp y3temp vectortemp r0temp
            
            x0temp=num2str(xyAlreadyCell{kpart}(3,1));
            y0temp=num2str(xyAlreadyCell{kpart}(3,2));
            r0temp=xyAlreadyCell{kpart}(3,3);
            vectortemp=[xBound;yBound]-[xyAlreadyCell{kpart}(3,1);xyAlreadyCell{kpart}(3,2)];
            vectortemp=(dot(vectortemp,[1;1])/abs(dot(vectortemp,[1;1]))*vectortemp/norm(vectortemp));
            x1temp=num2str(xyAlreadyCell{kpart}(3,1)+r0temp*vectortemp(1));
            y1temp=num2str(xyAlreadyCell{kpart}(3,2)+r0temp*vectortemp(2));
            fprintf(fID,'%s \n',['s1.ArcByCenterEnds(center=(',x0temp,', ',y0temp,'), point1=(',x1temp,', ',y1temp,'), point2=(',x1temp,', ',y1temp,'), direction=COUNTERCLOCKWISE)']);
            clear x0temp y0temp x1temp y1temp vectortemp 
            x0temp=xyAlreadyCell{kpart}(3,1);
            y0temp=xyAlreadyCell{kpart}(3,2);
            x2temp=CellEdgePoint{kpart}(4,1);
            y2temp=CellEdgePoint{kpart}(4,2);
            x3temp=CellEdgePoint{kpart}(1,1)+2*xBound;
            y3temp=CellEdgePoint{kpart}(1,2);
            vectortemp=0.5*[x2temp-x0temp;y2temp-y0temp]+0.5*[x3temp-x0temp;y3temp-y0temp];
            vectortemp=vectortemp/norm(vectortemp);
            InterfacePointMatrix=[InterfacePointMatrix;[x0temp+vectortemp(1)*r0temp,y0temp+vectortemp(2)*r0temp]];
            clear x0temp y0temp x2temp y2temp x3temp y3temp vectortemp r0temp
            
            x0temp=num2str(xyAlreadyCell{kpart}(4,1));
            y0temp=num2str(xyAlreadyCell{kpart}(4,2));
            r0temp=xyAlreadyCell{kpart}(4,3);
            vectortemp=[-xBound;yBound]-[xyAlreadyCell{kpart}(4,1);xyAlreadyCell{kpart}(4,2)];
            vectortemp=(dot(vectortemp,[-1;1])/abs(dot(vectortemp,[-1;1]))*vectortemp/norm(vectortemp));
            x1temp=num2str(xyAlreadyCell{kpart}(4,1)+r0temp*vectortemp(1));
            y1temp=num2str(xyAlreadyCell{kpart}(4,2)+r0temp*vectortemp(2));
            fprintf(fID,'%s \n',['s1.ArcByCenterEnds(center=(',x0temp,', ',y0temp,'), point1=(',x1temp,', ',y1temp,'), point2=(',x1temp,', ',y1temp,'), direction=COUNTERCLOCKWISE)']);         
             clear x0temp y0temp x1temp y1temp vectortemp 
            x0temp=xyAlreadyCell{kpart}(4,1);
            y0temp=xyAlreadyCell{kpart}(4,2);
            x2temp=CellEdgePoint{kpart}(1,1);
            y2temp=CellEdgePoint{kpart}(1,2);
            x3temp=CellEdgePoint{kpart}(3,1);
            y3temp=CellEdgePoint{kpart}(3,2);
            vectortemp=0.5*[x2temp-x0temp;y2temp-y0temp]+0.5*[x3temp-x0temp;y3temp-y0temp];
            vectortemp=vectortemp/norm(vectortemp);
            InterfacePointMatrix=[InterfacePointMatrix;[x0temp+vectortemp(1)*r0temp,y0temp+vectortemp(2)*r0temp]];
            clear x0temp y0temp x2temp y2temp x3temp y3temp vectortemp r0temp
            
        elseif ismember(kpart,xCrossArray) % x edge fiber
            x0temp=num2str(xyAlreadyCell{kpart}(1,1));
            y0temp=num2str(xyAlreadyCell{kpart}(1,2));
            r0temp=xyAlreadyCell{kpart}(1,3);
            x1temp=num2str(xyAlreadyCell{kpart}(1,1)-r0temp);
            y1temp=num2str(xyAlreadyCell{kpart}(1,2));
            fprintf(fID,'%s \n',['s1.ArcByCenterEnds(center=(',x0temp,', ',y0temp,'), point1=(',x1temp,', ',y1temp,'), point2=(',x1temp,', ',y1temp,'), direction=CLOCKWISE)']);
            InterfacePointMatrix=[InterfacePointMatrix;[xyAlreadyCell{kpart}(1,1)+r0temp,xyAlreadyCell{kpart}(1,2)]];
                
            x0temp=num2str(xyAlreadyCell{kpart}(2,1));
            y0temp=num2str(xyAlreadyCell{kpart}(2,2));
            r0temp=xyAlreadyCell{kpart}(2,3);
            x1temp=num2str(xyAlreadyCell{kpart}(1,1)+2*xBound+r0temp);
            y1temp=num2str(xyAlreadyCell{kpart}(1,2));
            fprintf(fID,'%s \n',['s1.ArcByCenterEnds(center=(',x0temp,', ',y0temp,'), point1=(',x1temp,', ',y1temp,'), point2=(',x1temp,', ',y1temp,'), direction=COUNTERCLOCKWISE)']);    
            InterfacePointMatrix=[InterfacePointMatrix;[xyAlreadyCell{kpart}(1,1)+2*xBound-r0temp,xyAlreadyCell{kpart}(1,2)]];
        elseif ismember(kpart,yCrossArray) % y edge fiber
            x0temp=num2str(xyAlreadyCell{kpart}(2,1));
            y0temp=num2str(xyAlreadyCell{kpart}(2,2));
            r0temp=xyAlreadyCell{kpart}(2,3);
            x1temp=num2str(xyAlreadyCell{kpart}(2,1));
            y1temp=num2str(xyAlreadyCell{kpart}(2,2)+r0temp);
            fprintf(fID,'%s \n',['s1.ArcByCenterEnds(center=(',x0temp,', ',y0temp,'), point1=(',x1temp,', ',y1temp,'), point2=(',x1temp,', ',y1temp,'), direction=COUNTERCLOCKWISE)']);
            InterfacePointMatrix=[InterfacePointMatrix;[xyAlreadyCell{kpart}(2,1),xyAlreadyCell{kpart}(2,2)-r0temp]];
            x0temp=num2str(xyAlreadyCell{kpart}(1,1));
            y0temp=num2str(xyAlreadyCell{kpart}(1,2));
            r0temp=xyAlreadyCell{kpart}(1,3);
            x1temp=num2str(xyAlreadyCell{kpart}(1,1));
            y1temp=num2str(xyAlreadyCell{kpart}(1,2)-r0temp);
            fprintf(fID,'%s \n',['s1.ArcByCenterEnds(center=(',x0temp,', ',y0temp,'), point1=(',x1temp,', ',y1temp,'), point2=(',x1temp,', ',y1temp,'), direction=CLOCKWISE)']);    
            InterfacePointMatrix=[InterfacePointMatrix;[xyAlreadyCell{kpart}(1,1),xyAlreadyCell{kpart}(1,2)+r0temp]];
        end
    end
    
    % Finish Plot figure1
    plot(InterfacePointMatrix(:,1),InterfacePointMatrix(:,2),'LineStyle','none','Marker','*','MarkerEdgeColor','b'); hold on;
    hold off;
    grid on;
    axis equal
    saveas(gcf,[dirname,'\model',num2str(ikpj),'.fig']);
    saveas(gcf,[dirname,'\model',num2str(ikpj),'.jpg']);
    
    
    fprintf(fID,'%s \n',['p = mdb.models[''Model-1''].parts[''Part-1'']']);
    fprintf(fID,'%s \n',['f = p.faces']);
    fprintf(fID,'%s \n',['pickedFaces = f.findAt(((0,0,0),),)']);
    fprintf(fID,'%s \n',['e1, d2 = p.edges, p.datums']);
    fprintf(fID,'%s \n',['p.PartitionFaceBySketch(faces=pickedFaces, sketch=s1)']);
    fprintf(fID,'%s \n',['s1.unsetPrimaryObject()']);
    fprintf(fID,'%s \n',['del mdb.models[''Model-1''].sketches[''__profile__'']']);

    % Material Properties
    fprintf(fID,'%s \n',['# Material Properties']);
    fprintf(fID,'%s \n',['session.viewports[''Viewport: 1''].setValues(displayedObject=None)']);
    fprintf(fID,'%s \n',['session.viewports[''Viewport: 1''].partDisplay.setValues(sectionAssignments=ON,engineeringFeatures=ON)']);
    fprintf(fID,'%s \n',['session.viewports[''Viewport: 1''].partDisplay.geometryOptions.setValues(referenceRepresentation=OFF)']);
    % material Fiber & Matrix
    fprintf(fID,'%s \n',['# Fiber']);
    fprintf(fID,'%s \n',['mdb.models[''Model-1''].Material(name=''Fiber'')']);
    fprintf(fID,'%s \n',['mdb.models[''Model-1''].materials[''Fiber''].Elastic(table=((230000.0, 0.3), ))']);
    fprintf(fID,'%s \n',['# Matrix']);
    fprintf(fID,'%s \n',['mdb.models[''Model-1''].Material(name=''Matrix'')']);
    fprintf(fID,'%s \n',['mdb.models[''Model-1''].materials[''Matrix''].Elastic(table=((3760.0, 0.39), ))']);
    % create sections
    fprintf(fID,'%s \n',['# Create Sections']);
    for kpart=1:nPart
        fprintf(fID,'%s \n',['mdb.models[''Model-1''].HomogeneousSolidSection(name=''Fiber',num2str(kpart),''', material=''Fiber'',thickness=None)']);
    end
    fprintf(fID,'%s \n',['mdb.models[''Model-1''].HomogeneousSolidSection(name=''Matrix'', material=''Matrix'',thickness=None)']);
    % Assign Sections
    fprintf(fID,'%s \n',['# Assign Sections']);
    fprintf(fID,'%s \n',['p = mdb.models[''Model-1''].parts[''Part-1'']']);
    fprintf(fID,'%s \n',['f = p.faces']);
    MaterialPointCell=cell(nPart+1,1);
    for kpart=1:nPart
        fprintf(fID,'%s \n',['# Fiber Section ',num2str(kpart)]);
        if ~ismember(kpart,[cornerCrossArray,xCrossArray,yCrossArray])  % inner fiber
            x0temp=xyAlreadyCell{kpart}(1,1);
            y0temp=xyAlreadyCell{kpart}(1,2);
            fprintf(fID,'%s \n',['faces=f.findAt((( ',num2str(x0temp),', ',num2str(y0temp),',0.0),),)']);
            MaterialPointCell{kpart}=[x0temp,y0temp];
            clear x0temp y0temp
        elseif ismember(kpart,cornerCrossArray) % Corner fiber
            referencePoint=zeros(4,2);
            x0temp=CellEdgePoint{kpart}(2,1);
            y0temp=CellEdgePoint{kpart}(2,2);
            x1temp=CellEdgePoint{kpart}(3,1);
            y1temp=CellEdgePoint{kpart}(3,2)-2*yBound;
            x2temp=0.5*x0temp+0.5*x1temp;
            y2temp=0.5*y0temp+0.5*y1temp;
            referencePoint(1,:)=[x2temp,y2temp];
            clear x0temp y0temp x1temp y1temp x2temp y2temp
            x0temp=CellEdgePoint{kpart}(4,1);
            y0temp=CellEdgePoint{kpart}(4,2)-2*yBound;
            x1temp=CellEdgePoint{kpart}(2,1)+2*xBound;
            y1temp=CellEdgePoint{kpart}(2,2);
            x2temp=0.5*x0temp+0.5*x1temp;
            y2temp=0.5*y0temp+0.5*y1temp;
            referencePoint(2,:)=[x2temp,y2temp];
            clear x0temp y0temp x1temp y1temp x2temp y2temp vectortemp
            x0temp=CellEdgePoint{kpart}(4,1);
            y0temp=CellEdgePoint{kpart}(4,2);
            x1temp=CellEdgePoint{kpart}(1,1)+2*xBound;
            y1temp=CellEdgePoint{kpart}(1,2);
            x2temp=0.5*x0temp+0.5*x1temp;
            y2temp=0.5*y0temp+0.5*y1temp;
            referencePoint(3,:)=[x2temp,y2temp];
            clear x0temp y0temp x1temp y1temp x2temp y2temp vectortemp
            x0temp=CellEdgePoint{kpart}(1,1);
            y0temp=CellEdgePoint{kpart}(1,2);
            x1temp=CellEdgePoint{kpart}(3,1);
            y1temp=CellEdgePoint{kpart}(3,2);
            x2temp=0.5*x0temp+0.5*x1temp;
            y2temp=0.5*y0temp+0.5*y1temp;
            referencePoint(4,:)=[x2temp,y2temp];
            clear x0temp y0temp x1temp y1temp x2temp y2temp vectortemp
            MaterialPointCell{kpart}=referencePoint;
            x1temp=num2str(referencePoint(1,1));    y1temp=num2str(referencePoint(1,2));
            x2temp=num2str(referencePoint(2,1));    y2temp=num2str(referencePoint(2,2));
            x3temp=num2str(referencePoint(3,1));    y3temp=num2str(referencePoint(3,2));
            x4temp=num2str(referencePoint(4,1));    y4temp=num2str(referencePoint(4,2));
            fprintf(fID,'%s \n',['faces=f.findAt(((',x1temp,',',y1temp,',0.0),),((',x2temp,',',y2temp,',0.0),),((',x3temp,',',y3temp,',0.0),),((',x4temp,',',y4temp,',0.0),),)']);
            
            clear x1temp x2temp x3temp x4temp y1temp y2temp y3temp y4temp referencePoint
        elseif ismember(kpart,xCrossArray) % x edge fiber
            referencePoint=zeros(2,2);
            x0temp=xyAlreadyCell{kpart}(1,1);
            y0temp=xyAlreadyCell{kpart}(1,2);
            r0temp=xyAlreadyCell{kpart}(1,3);
            x1temp=x0temp+r0temp;
            y1temp=y0temp;
            x2temp=0.5*x1temp+0.5*(-xBound);
            y2temp=y1temp;
            referencePoint(1,:)=[x2temp,y2temp];
            clear x0temp y0temp x1temp y1temp x2temp y2temp 
            x0temp=xyAlreadyCell{kpart}(2,1);
            y0temp=xyAlreadyCell{kpart}(2,2);
            r0temp=xyAlreadyCell{kpart}(2,3);
            x1temp=x0temp-r0temp;
            y1temp=y0temp;
            x2temp=0.5*x1temp+0.5*xBound;
            y2temp=y1temp;
            referencePoint(2,:)=[x2temp,y2temp];
            clear x0temp y0temp x1temp y1temp x2temp y2temp 
            MaterialPointCell{kpart}=referencePoint;
            x1temp=num2str(referencePoint(1,1));    y1temp=num2str(referencePoint(1,2));
            x2temp=num2str(referencePoint(2,1));    y2temp=num2str(referencePoint(2,2));
            fprintf(fID,'%s \n',['faces=f.findAt(((',x1temp,',',y1temp,',0.0),),((',x2temp,',',y2temp,',0.0),),)']);
            clear x1temp x2temp y1temp y2temp referencePoint
        elseif ismember(kpart,yCrossArray) % y edge fiber
            referencePoint=zeros(2,2);
            x0temp=xyAlreadyCell{kpart}(1,1);
            y0temp=xyAlreadyCell{kpart}(1,2);
            r0temp=xyAlreadyCell{kpart}(1,3);
            x1temp=x0temp;
            y1temp=y0temp+r0temp;
            x2temp=x1temp;
            y2temp=0.5*y1temp+0.5*(-yBound);
            referencePoint(1,:)=[x2temp,y2temp];
            clear x0temp y0temp x1temp y1temp x2temp y2temp 
            x0temp=xyAlreadyCell{kpart}(2,1);
            y0temp=xyAlreadyCell{kpart}(2,2);
            r0temp=xyAlreadyCell{kpart}(2,3);
            x1temp=x0temp;
            y1temp=y0temp-r0temp;
            x2temp=x1temp;
            y2temp=0.5*y1temp+0.5*yBound;
            referencePoint(2,:)=[x2temp,y2temp];
            clear x0temp y0temp x1temp y1temp x2temp y2temp 
            MaterialPointCell{kpart}=referencePoint;
            x1temp=num2str(referencePoint(1,1));    y1temp=num2str(referencePoint(1,2));
            x2temp=num2str(referencePoint(2,1));    y2temp=num2str(referencePoint(2,2));
            fprintf(fID,'%s \n',['faces=f.findAt(((',x1temp,',',y1temp,',0.0),),((',x2temp,',',y2temp,',0.0),),)']);
            clear x1temp x2temp y1temp y2temp referencePoint
        end
        fprintf(fID,'%s \n',['region = p.Set(faces=faces, name=''Fiber',num2str(kpart),''')']);
        fprintf(fID,'%s \n',['p = mdb.models[''Model-1''].parts[''Part-1'']']);
        fprintf(fID,'%s \n',['p.SectionAssignment(region=region, sectionName=''Fiber',num2str(kpart),''', offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='''', thicknessAssignment=FROM_SECTION)']);
    end
    
    % Assign Matrix Section
    % Get a Point in Matrix  
    x1temp=xyAlreadyCell{2}(1,1);
    y1temp=xyAlreadyCell{2}(1,2);
    r1temp=xyAlreadyCell{2}(1,3);
    middleIter=0;
    xyAlreadytemp=[];
    for kpart=1:nPart
        leng1=size(xyAlreadyCell{kpart},1);
        xyAlreadytemp=[xyAlreadytemp;[kpart*ones(leng1,1),xyAlreadyCell{kpart}]];
        clear leng1
    end
    ifSucess=false;
    for kpart=1:nPart
        x0temp=xyAlreadyCell{kpart}(1,1);
        y0temp=xyAlreadyCell{kpart}(1,2);
        r0temp=xyAlreadyCell{kpart}(1,3);
        while true
            middleIter=middleIter+1;
            if middleIter>10000
                error('cannot find a matrix material point')
            end
            xmiddletemp=x0temp*r1temp/(r0temp+r1temp)+x1temp*r0temp/(r0temp+r1temp);
            ymiddletemp=y0temp*r1temp/(r1temp+r1temp)+y1temp*r0temp/(r0temp+r1temp);
            
            NormalizeDistanceArray=sqrt((xmiddletemp-xyAlreadytemp(:,2)).^2+(ymiddletemp-xyAlreadytemp(:,3)).^2)./xyAlreadytemp(:,4);
            InPart=find(NormalizeDistanceArray<1.0);
            if isempty(InPart)&&(abs(xmiddletemp)<xBound)&&(abs(ymiddletemp)<yBound)
                Matrixpoint=[xmiddletemp,ymiddletemp];
                clear x0temp y0temp r0temp x1temp y1temp r1temp NormalizeDistanceArray InPart xmiddletemp ymiddletemp xyAlreadytemp
                ifSucess=true;    
                break;
            else
                x1temp=xyAlreadyCell{xyAlreadytemp(InPart(1),1)}(1,1);
                y1temp=xyAlreadyCell{xyAlreadytemp(InPart(1),1)}(1,2);
                r1temp=xyAlreadyCell{xyAlreadytemp(InPart(1),1)}(1,3);
            end
        end
        if ifSucess
            break; 
        end
    end
    fprintf(fID,'%s \n',['# Matrix Section']);
    fprintf(fID,'%s \n',['faces=f.findAt((( ',num2str(Matrixpoint(1)),', ',num2str(Matrixpoint(2)),',0.0),),)']);
    MaterialPointCell{nPart+1}=[Matrixpoint(1),Matrixpoint(2)];
    fprintf(fID,'%s \n',['region = p.Set(faces=faces, name=''Matrix'')']);
    fprintf(fID,'%s \n',['p = mdb.models[''Model-1''].parts[''Part-1'']']);
    fprintf(fID,'%s \n',['p.SectionAssignment(region=region, sectionName=''Matrix'', offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='''', thicknessAssignment=FROM_SECTION)']);
    % Mesh the Part
    fprintf(fID,'%s \n',['session.viewports[''Viewport: 1''].setValues(displayedObject=None)']);
    fprintf(fID,'%s \n',['p = mdb.models[''Model-1''].parts[''Part-1'']']);
    fprintf(fID,'%s \n',['session.viewports[''Viewport: 1''].setValues(displayedObject=p)']);
    fprintf(fID,'%s \n',['session.viewports[''Viewport: 1''].partDisplay.setValues(sectionAssignments=OFF,engineeringFeatures=OFF, mesh=ON)']);
    fprintf(fID,'%s \n',['session.viewports[''Viewport: 1''].partDisplay.meshOptions.setValues(meshTechnique=ON)']);
    fprintf(fID,'%s \n',['mdb.models[''Model-1''].parts[''Part-1''].setValues(geometryRefinement=EXTRA_FINE)']);
    fprintf(fID,'%s \n',['p = mdb.models[''Model-1''].parts[''Part-1'']']);
    fprintf(fID,'%s \n',['e = p.edges']);
    % Seed the Edges
    % Inner
    fprintf(fID,'%s \n',['# Inner Interface Seed ']);
    leng1=size(InterfacePointMatrix,1);
    for kedge=1:leng1
        fprintf(fID,'%s \n',['PickedgeDic=e.getClosest(coordinates=((',num2str(InterfacePointMatrix(kedge,1)),', ',num2str(InterfacePointMatrix(kedge,2)),', 0.0),),)']);
        fprintf(fID,'%s \n',['pickedEdges=(PickedgeDic[0][0],)']);
        fprintf(fID,'%s \n',['p.seedEdgeBySize(edges=pickedEdges, size=',num2str(lengElem0),', deviationFactor=0.1,constraint=FINER)']);
    end
    % Left
    fprintf(fID,'%s \n',['# Left Boundary Seed ']);
    leng1=length(EdgeXmiddle);
    for kedge=1:leng1
        fprintf(fID,'%s \n',['pickedEdges=e.findAt(((',num2str(-xBound),',',num2str(EdgeXmiddle(kedge)),',0.0),),)']);
        fprintf(fID,'%s \n',['p.seedEdgeBySize(edges=pickedEdges, size=',num2str(lengElem0),', deviationFactor=0.1,constraint=FINER)']);
    end
    % Right
    fprintf(fID,'%s \n',['# Right Boundary Seed ']);
    for kedge=1:leng1
        fprintf(fID,'%s \n',['pickedEdges=e.findAt(((',num2str(xBound),',',num2str(EdgeXmiddle(kedge)),',0.0),),)']);
        fprintf(fID,'%s \n',['p.seedEdgeBySize(edges=pickedEdges, size=',num2str(lengElem0),', deviationFactor=0.1,constraint=FINER)']);
    end
    % Up
    fprintf(fID,'%s \n',['# Up Boundary Seed ']);
    leng1=length(EdgeYmiddle);
    for kedge=1:leng1
        fprintf(fID,'%s \n',['pickedEdges=e.findAt(((',num2str(EdgeYmiddle(kedge)),',',num2str(yBound),',0.0),),)']);
        fprintf(fID,'%s \n',['p.seedEdgeBySize(edges=pickedEdges, size=',num2str(lengElem0),', deviationFactor=0.1,constraint=FINER)']);
    end
    % Down
    fprintf(fID,'%s \n',['# Down Boundary Seed ']);
    for kedge=1:leng1
        fprintf(fID,'%s \n',['pickedEdges=e.findAt(((',num2str(EdgeYmiddle(kedge)),',',num2str(-yBound),',0.0),),)']);
        fprintf(fID,'%s \n',['p.seedEdgeBySize(edges=pickedEdges, size=',num2str(lengElem0),', deviationFactor=0.1,constraint=FINER)']);
    end
    fprintf(fID,'%s \n',['p.seedEdgeBySize(edges=pickedEdges, size=',num2str(lengElem0),', deviationFactor=0.1,constraint=FINER)  ']);
    fprintf(fID,'%s \n',['p.seedPart(size=',num2str(lengElem1),', deviationFactor=0.1, minSizeFactor=0.1) ']);
    fprintf(fID,'%s \n',['p = mdb.models[''Model-1''].parts[''Part-1'']']);
    fprintf(fID,'%s \n',['f = p.faces ']);
    for kpart=1:nPart+1
        leng1=size(MaterialPointCell{kpart},1);
        for ksub=1:leng1
            fprintf(fID,'%s \n',['pickedRegions = f.findAt(((',num2str(MaterialPointCell{kpart}(ksub,1)),',',num2str(MaterialPointCell{kpart}(ksub,2)),',0.),),)']);
            fprintf(fID,'%s \n',['p.setMeshControls(regions=pickedRegions, elemShape=QUAD)']);
        end
    end
    fprintf(fID,'%s \n',['elemType1 = mesh.ElemType(elemCode=CPS4, elemLibrary=STANDARD)']);
    fprintf(fID,'%s \n',['elemType2 = mesh.ElemType(elemCode=CPS3, elemLibrary=STANDARD)']);
    fprintf(fID,'%s \n',['p = mdb.models[''Model-1''].parts[''Part-1'']']);
    fprintf(fID,'%s \n',['f = p.faces']);
    for kpart=1:nPart+1
        leng1=size(MaterialPointCell{kpart},1);
        for ksub=1:leng1
            fprintf(fID,'%s \n',['faces = f.findAt(((',num2str(MaterialPointCell{kpart}(ksub,1)),',',num2str(MaterialPointCell{kpart}(ksub,2)),',0.),),)']);
            fprintf(fID,'%s \n',['pickedRegions =(faces, )']);
            fprintf(fID,'%s \n',['p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))']);
        end
    end
    fprintf(fID,'%s \n',['p = mdb.models[''Model-1''].parts[''Part-1'']']);
    fprintf(fID,'%s \n',['p.generateMesh()']);
    fprintf(fID,'%s \n',['cmap=session.viewports[''Viewport: 1''].colorMappings[''Section'']']);
    fprintf(fID,'%s \n',['session.viewports[''Viewport: 1''].setColor(colorMapping=cmap)']);
    fprintf(fID,'%s \n',['session.viewports[''Viewport: 1''].disableMultipleColors()']);
    fprintf(fID,'%s \n',['session.printToFile(fileName=''Model',num2str(ikpj),''',format=PNG, canvasObjects=(session.viewports[''Viewport: 1''], ))']);
    % Assembly
    fprintf(fID,'%s \n',['a = mdb.models[''Model-1''].rootAssembly']);
    fprintf(fID,'%s \n',['a.DatumCsysByDefault(CARTESIAN)']);
    fprintf(fID,'%s \n',['p = mdb.models[''Model-1''].parts[''Part-1'']']);
    fprintf(fID,'%s \n',['a.Instance(name=''Part-1-1'', part=p, dependent=ON)']);
    
    % .inp writer
    fprintf(fID,'%s \n',['mdb.Job(name=''',['Model',num2str(ikpj)],''', model=''Model-1'')']);
    fprintf(fID,'%s \n',['mdb.jobs[''',['Model',num2str(ikpj)],'''].writeInput(consistencyChecking=OFF)']);
    fprintf(fID,'%s \n',['mdb.saveAs(pathName=''',fileName,''')']);
    
    fclose(fID);
    % 
    % .py Run
    cd(dirname);
    system([' abaqus cae noGUI=',fileName,'.py']);
    cd ..;

    %% Read Old .inp File
    fID=fopen([dirname,'\Model',num2str(ikpj),'.inp'],'r');
    if fID<0
        continue;
    end
    StrInfor=textscan(fID, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s','Delimiter',',');
    fclose(fID);
    Str1Cell=StrInfor{1};
    Str3Cell=StrInfor{3};
    clear StrInfor
    nodelineLoc=find(strcmp(Str1Cell,'*Node'));     % Colume of the Start of Nodes
    elemlineLoc=find(strcmp(Str1Cell,'*Element'));  % Colume of the Start of Elements
    elsetlineLoc=find(strcmp(Str1Cell,'*Elset'));   % Colume of the Start of Elsets
    elsetGeneratelineLoc=find(strcmp(Str3Cell,'generate')); % Colume of the Elset defined by generate method
    clear Str1Cell Str3Cell
    
    
    leng1=length(elsetGeneratelineLoc);
    generateArray=zeros(1,leng1);
    for kset=1:leng1
        generateArray(kset)=find(elsetlineLoc==elsetGeneratelineLoc(kset));
    end
    
    % Get Node Information
    fID=fopen([dirname,'\Model',num2str(ikpj),'.inp'],'r');
    StrInfor=textscan(fID, '%f %f %f ','Delimiter',',','HeaderLines',nodelineLoc);
    fclose(fID);
    nodeInfor=[StrInfor{1},StrInfor{2},StrInfor{3}];
    numNode=size(nodeInfor,1);  % number of total nodes
    numNodeInput=numNode;       % initial number
    clear StrInfor
    
    % Get Element Information
    fID=fopen([dirname,'\Model',num2str(ikpj),'.inp'],'r');
    StrInfor=textscan(fID, '%f %f %f %f %f','Delimiter',',','HeaderLines',elemlineLoc);
    fclose(fID);
    elemInfor=[StrInfor{1},StrInfor{2},StrInfor{3},StrInfor{4},StrInfor{5}];
    numElem=size(elemInfor,1);  % number of total elements
    numElemInput=numElem;       % initial number
    clear StrInfor
    elemSetInfor=cell(nPart+1,1);
    for kpart=1:nPart+1
        if ismember(kpart,generateArray)
            fID=fopen([dirname,'\Model',num2str(ikpj),'.inp'],'r');
            StrInfor=textscan(fID, '%f %f %f','Delimiter',',','HeaderLines',elsetlineLoc(kpart));
            fclose(fID);
            elemSetInfor{kpart}=transpose(StrInfor{1}:StrInfor{3}:StrInfor{2});
            clear StrInfor
        else            
            fID=fopen([dirname,'\Model',num2str(ikpj),'.inp'],'r');
            StrInfor=textscan(fID, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','Delimiter',',','HeaderLines',elsetlineLoc(kpart));
            fclose(fID);
            elemSetInforTemp=[StrInfor{1};StrInfor{2};StrInfor{3};StrInfor{4};StrInfor{5};StrInfor{6};StrInfor{7};StrInfor{8};StrInfor{9};StrInfor{10};StrInfor{11};StrInfor{12};StrInfor{13};StrInfor{14};StrInfor{15};StrInfor{16}];
            elemSetInforTemp(isnan(elemSetInforTemp))=[];
            elemSetInfor{kpart}=sort(elemSetInforTemp);
            clear StrInfor elemSetInforTemp
        end
    end
    elemInforCell=cell(nPart,1);
    for kpart=1:nPart+1
        elemInforCell{kpart}=elemInfor(elemSetInfor{kpart},:);
    end
    
    % find all outer edge
    % The Local Nodal ID order
    %   4 ------------- 3
    %   |               |
    %   |    element    |
    %   |     CPS4      |
    %   |               |
    %   1 ------------- 2
    FaceInforMatrix=[[elemInfor(:,2),elemInfor(:,3)];[elemInfor(:,3),elemInfor(:,4)];[elemInfor(:,4),elemInfor(:,5)];[elemInfor(:,5),elemInfor(:,2)]];
    FaceInforMatrix=sort(FaceInforMatrix,2); % sort each row
    FaceInforMatrixUnique=unique(FaceInforMatrix,'rows');   % unique the faces
    [FaceTemp,ia,ib]=intersect(FaceInforMatrix,FaceInforMatrixUnique,'rows');
    FaceInforMatrixInner=FaceInforMatrix;
    FaceInforMatrixInner(ia,:)=[];
    FaceInforMatrixOuter = setdiff(FaceInforMatrixUnique,FaceInforMatrixInner,'rows');
    
    FreeNodeInforMatrix=[nodeInfor(FaceInforMatrixOuter(:,1),2:end),nodeInfor(FaceInforMatrixOuter(:,2),2:end)];                                    % 1-4, record the coord of 2 point
    FreeNodeInforMatrix=[FreeNodeInforMatrix,transpose(cross(transpose([FreeNodeInforMatrix(:,1:2)-FreeNodeInforMatrix(:,3:4),zeros(size(FreeNodeInforMatrix,1),1)]),kron(ones(1,size(FreeNodeInforMatrix,1),1),[0;0;1])))];    % 5-7, calculate 1 normal vector
    FreeNodeInforMatrix(:,end)=[];      % Set z axis free
    FreeNodeInforMatrix=[FreeNodeInforMatrix,sqrt(FreeNodeInforMatrix(:,5).^2+FreeNodeInforMatrix(:,6).^2)];                 % 7, the length of the normal vector
    FreeNodeInforMatrix(:,5)=FreeNodeInforMatrix(:,5)./FreeNodeInforMatrix(:,7);
    FreeNodeInforMatrix(:,6)=FreeNodeInforMatrix(:,6)./FreeNodeInforMatrix(:,7);
    
    FreeNodeInforMatrix=[FreeNodeInforMatrix,FreeNodeInforMatrix(:,1).*FreeNodeInforMatrix(:,5)+FreeNodeInforMatrix(:,2).*FreeNodeInforMatrix(:,6)];  % 8, if the noraml vector point to out
    FreeNodeInforMatrix(:,8)= FreeNodeInforMatrix(:,8)./abs(FreeNodeInforMatrix(:,8));
    FreeNodeInforMatrix(:,5)=FreeNodeInforMatrix(:,5).*FreeNodeInforMatrix(:,8);                                                                     % make the normal from inner to outer
    FreeNodeInforMatrix(:,6)=FreeNodeInforMatrix(:,6).*FreeNodeInforMatrix(:,8);

    FreeNodeInforMatrix=[FreeNodeInforMatrix,round(2*FreeNodeInforMatrix(:,5)+FreeNodeInforMatrix(:,6))];                             % 9, the symbol of the different face
    FreeNodeInforMatrix=[FreeNodeInforMatrix,nodeInfor(FaceInforMatrixOuter(:,1),1),nodeInfor(FaceInforMatrixOuter(:,2),1)];          % 10-11,ID  of nodes
    
    FaceXP=find(FreeNodeInforMatrix(:,9)==2);
    FaceXN=find(FreeNodeInforMatrix(:,9)==-2);
    FaceYP=find(FreeNodeInforMatrix(:,9)==1);
    FaceYN=find(FreeNodeInforMatrix(:,9)==-1);
     
    NodeXP=[FreeNodeInforMatrix(FaceXP,10);FreeNodeInforMatrix(FaceXP,11)];
    NodeXP=unique(NodeXP);
    NodeXPinfor=nodeInfor(NodeXP,:);
    NodeXPinfor=sortrows(NodeXPinfor,3);
    NodeXP=NodeXPinfor(:,1);
    clear NodeXPinfor
    NodeXN=[FreeNodeInforMatrix(FaceXN,10);FreeNodeInforMatrix(FaceXN,11)];
    NodeXN=unique(NodeXN);
    NodeXNinfor=nodeInfor(NodeXN,:);
    NodeXNinfor=sortrows(NodeXNinfor,3);
    NodeXN=NodeXNinfor(:,1);
    clear NodeXNinfor
    NodeYP=[FreeNodeInforMatrix(FaceYP,10);FreeNodeInforMatrix(FaceYP,11)];
    NodeYP=unique(NodeYP);
    NodeYPinfor=nodeInfor(NodeYP,:);
    NodeYPinfor=sortrows(NodeYPinfor,2);
    NodeYP=NodeYPinfor(:,1);
    clear NodeYPinfor
    NodeYN=[FreeNodeInforMatrix(FaceYN,10);FreeNodeInforMatrix(FaceYN,11)];
    NodeYN=unique(NodeYN);
    NodeYNinfor=nodeInfor(NodeYN,:);
    NodeYNinfor=sortrows(NodeYNinfor,2);
    NodeYN=NodeYNinfor(:,1);
    clear NodeYNinfor
    
    clear FreeNodeInforMatrix FaceXP FaceXN FaceYP FaceYN
    
    if (size(NodeXP,1)~=size(NodeXN,1))||(size(NodeYP,1)~=size(NodeYN,1))
        warning(['This seed(',num2str(ikpj),') is not suitable to generate model: Abaqus mesh error']);
        system(['del ',dirname,'\Part2D_',num2str(ikpj),'.cae']);
        system(['del ',dirname,'\Part2D_',num2str(ikpj),'.jnl']);
        system(['del ',dirname,'\Part2D_',num2str(ikpj),'.py']);
        system(['del ',dirname,'\model',num2str(ikpj),'.fig']);
        system(['del ',dirname,'\model',num2str(ikpj),'.jpg']);
        system(['del ',dirname,'\Model',num2str(ikpj),'.inp']);
        system(['del ',dirname,'\Model',num2str(ikpj),'.png']);
        continue;
    end
    
    % Get the nodes at boundary btw Fiber & Matrix
    
    
    if addInterface||YLmemb
        elemInforMatrix_Part=[];    % elements in particles/fibers
        elemInforMatrix_Matrix=[];  % elements in the matrix
        for kpart = 1:nPart
            elemInforMatrix_Part=[elemInforMatrix_Part;elemInforCell{kpart}];
        end
        elemInforMatrix_Matrix=elemInforCell{1+nPart};
        Face_Part=[[elemInforMatrix_Part(:,2),elemInforMatrix_Part(:,3)];[elemInforMatrix_Part(:,3),elemInforMatrix_Part(:,4)];[elemInforMatrix_Part(:,4),elemInforMatrix_Part(:,5)];[elemInforMatrix_Part(:,5),elemInforMatrix_Part(:,2)]];
        Face_Part=sort(Face_Part,2);
        Face_Matrix=[[elemInforMatrix_Matrix(:,2),elemInforMatrix_Matrix(:,3)];[elemInforMatrix_Matrix(:,3),elemInforMatrix_Matrix(:,4)];[elemInforMatrix_Matrix(:,4),elemInforMatrix_Matrix(:,5)];[elemInforMatrix_Matrix(:,5),elemInforMatrix_Matrix(:,2)]];
        Face_Matrix=sort(Face_Matrix,2);
        
        FaceShare=intersect(Face_Part,Face_Matrix,'rows');  % Get the sharing Face in the mesh
        FaceShare_Num=size(FaceShare,1);                    % number of edges on interface
        FaceShare_Node=[FaceShare(:,1);FaceShare(:,2)];
        FaceShare_Node_Unique=unique(FaceShare_Node);       % Get all nodes lying on the interface
        FaceShare_Node_Num=size(FaceShare_Node_Unique,1);   % number of sharing nodes
        
    end 
    if YLmemb
        MemEle_Matrix=FaceShare;
        MemEle_Matrix=[transpose((numElem+1):(FaceShare_Num+numElem)),MemEle_Matrix];
        numElem=numElem+FaceShare_Num;
    end
    if addInterface
        NodeNew=transpose((numNode+1):(FaceShare_Node_Num+numNode));     % Give blank number for new nodes
        NodePair=[FaceShare_Node_Unique,NodeNew];                        % record new nodes
        nodeInfor(NodeNew,:)=[NodeNew,nodeInfor(FaceShare_Node_Unique,2:end)]; 
        NodePair_Old=NodePair(:,1);
        NodePair_New=NodePair(:,2);
        numNode=numNode+FaceShare_Node_Num;                              % Update the nodes total number
        
        FaceShare_Node1=FaceShare(:,1);
        FaceShare_Node2=FaceShare(:,2);
        
        CohEle_Matrix=NaN*ones(FaceShare_Num,2);    % the matrix record the new COHESIVE element
        for fNum=1:FaceShare_Num
            CohEle_Matrix(fNum,:)=[NodePair_New(find(NodePair_Old==FaceShare_Node2(fNum))),NodePair_New(find(NodePair_Old==FaceShare_Node1(fNum)))]; % the new nodes in cohesive element
        end
        CohEle_Matrix=[transpose((numElem+1):(FaceShare_Num+numElem)),FaceShare,CohEle_Matrix];
        numElem=numElem+FaceShare_Num;  % Update the elements total number
        
        elemInforMatrix_MatrixnoID=elemInforMatrix_Matrix(:,2:end); % replace the old nodes in matrix with new nodes
        for nNum=1:FaceShare_Node_Num    % Relpace the old nodes
            [ReplaceLoc1,ReplaceLoc2]=find(elemInforMatrix_MatrixnoID==NodePair(nNum,1));
            for kReplace=1:size(ReplaceLoc1,1)
                elemInforMatrix_Matrix(ReplaceLoc1(kReplace),ReplaceLoc2(kReplace)+1)=NodePair(nNum,2);
            end
        end
        elemInforCell{nPart+1}=elemInforMatrix_Matrix;
        
        NodeXPnew=NodeXP;
        NodeXNnew=NodeXN;
        NodeYPnew=NodeYP;
        NodeYNnew=NodeYN;
        % x
        [InterPhaseXNtemp,invx0,iNPO]=intersect(NodeXP,NodePair_Old);  % find the all 'Boundary Interface Nodes' and return their location in <nvx0>
        OrderNew=sortrows([iNPO,invx0],2);                % Sort the new nodes by the asccending order of old nodes
        NodeXPnew=[NodeXPnew;NodePair_New(OrderNew(:,1))];                    % The new xneg Bounary Node Set
        
        [InterPhaseXPtemp,invx1,iNPO]=intersect(NodeXN,NodePair_Old);  % find the all 'Boundary Interface Nodes' and return their location in <nvx0>
        OrderNew=sortrows([iNPO,invx1],2);                % Sort the new nodes by the asccending order of old nodes
        NodeXNnew=[NodeXNnew;NodePair_New(OrderNew(:,1))];                    % The new xpos Bounary Node Set
        
        % y
        [InterPhaseYNtemp,invy0,iNPO]=intersect(NodeYP,NodePair_Old);  % find the all 'Boundary Interface Nodes' and return their location in <nvx0>
        OrderNew=sortrows([iNPO,invy0],2);                % Sort the new nodes by the asccending order of old nodes
        NodeYPnew=[NodeYPnew;NodePair_New(OrderNew(:,1))];                    % The new xneg Bounary Node Set
        
        [InterPhaseYPtemp,invy1,iNPO]=intersect(NodeYN,NodePair_Old);  % find the all 'Boundary Interface Nodes' and return their location in <nvx0>
        OrderNew=sortrows([iNPO,invy1],2);                % Sort the new nodes by the asccending order of old nodes
        NodeYNnew=[NodeYNnew;NodePair_New(OrderNew(:,1))];                    % The new xpos Bounary Node Set
        
        NodeXP=NodeXPnew;
        NodeXN=NodeXNnew;
        NodeYP=NodeYPnew;
        NodeYN=NodeYNnew;
        
       
        clear NodeXPnew NodeXNnew NodeYPnew NodeYNnew OrderNew invy1 iNPO InterPhaseYPtemp
    end
    
    NodeInner=transpose(1:numNode);
    NodeInner=setdiff(NodeInner,NodeXP);
    NodeInner=setdiff(NodeInner,NodeXN);
    NodeInner=setdiff(NodeInner,NodeYP);
    NodeInner=setdiff(NodeInner,NodeYN);
    
    Node_A=intersect(NodeXN,NodeYN);
    Node_B=intersect(NodeXP,NodeYN);
    Node_C=intersect(NodeXP,NodeYP);
    Node_D=intersect(NodeXN,NodeYP);
    
%     NodeXN=setdiff(NodeXN,[Node_A;Node_D]);
    NodeXN(1)=[];
    NodeXN(end)=[];
%     NodeXP=setdiff(NodeXP,[Node_B;Node_C]);
    NodeXP(1)=[];
    NodeXP(end)=[];
%     NodeYN=setdiff(NodeYN,[Node_A;Node_B]);
    NodeYN(1)=[];
    NodeYN(end)=[];
%     NodeYP=setdiff(NodeYP,[Node_C;Node_D]);
    NodeYP(1)=[];
    NodeYP(end)=[];
    
    [~,Node_Alpha]=min((nodeInfor(NodeInner,2)-nodeInfor(Node_A,2)).^2+...
        (nodeInfor(NodeInner,3)-nodeInfor(Node_A,3)).^2);
    Node_Alpha=NodeInner(Node_Alpha,1);
    [~,Node_Beta]=min((nodeInfor(NodeInner,2)-nodeInfor(Node_B,2)).^2+...
        (nodeInfor(NodeInner,3)-nodeInfor(Node_B,3)).^2);
    Node_Beta=NodeInner(Node_Beta,1);
    
    
    %% WRITE NEW .INP FILE
    fileNameNew=[dirname,'\',fileName,'_PBC'];
    if YLmemb
        fileNameNew=[fileNameNew,'_YL'];
    end
    if addInterface
        fileNameNew=[fileNameNew,'_COH'];
    end
    fileNameNew=[fileNameNew,'.inp'];
    fID=fopen(fileNameNew,'w');
    %
    % heading
    %
    fprintf(fID,'%s\n',['*Heading']);
    abc=datestr(now);
    fprintf(fID,'%s\n',['ABAQUS job created by Peiliang Bian on ' abc(1:11) ' at ' abc(13:end)] );
    %
    % part
    %
    fprintf(fID,'%s\n',['**']);
    fprintf(fID,'%s\n',['** PART']);
    fprintf(fID,'%s\n',['**']);
    fprintf(fID,'%s\n',['*Part, name=Part-1']);
    fprintf(fID,'%s\n',['*Node']);
    fclose(fID);
    dlmwrite(fileNameNew,nodeInfor,'-append','precision',10);
    fID=fopen(fileNameNew,'a');
    fprintf(fID,'%s \n',['*NODE, nset=nodeXX']);
    fprintf(fID,'%s \n',[num2str(numNode+1) ', ,']);
    fprintf(fID,'%s \n',['*NODE, nset=nodeYY']);
    fprintf(fID,'%s \n',[num2str(numNode+2) ', ,']);
    fprintf(fID,'%s \n',['*NODE, nset=nodeXY']);
    fprintf(fID,'%s \n',[num2str(numNode+3) ', ,']);
    fclose(fID);
    for kpart=1:nPart
        fID=fopen(fileNameNew,'a');
        fprintf(fID,'%s\n',['*Element, type=CPE4, elset=PART-',num2str(kpart)]);
        fclose(fID);
        dlmwrite(fileNameNew,elemInforCell{kpart},'-append','precision',10);
    end
    fID=fopen(fileNameNew,'a');
    fprintf(fID,'%s\n',['*Element, type=CPE4, elset=MATRIX']);
    dlmwrite(fileNameNew,elemInforCell{end},'-append','precision',10);
    fclose(fID);
    if YLmemb
        fID=fopen(fileNameNew,'a');
        fprintf(fID,'%s\n',['*Element, type=T2D2, elset=YL']);
        dlmwrite(fileNameNew,MemEle_Matrix,'-append','precision',10);
        fclose(fID);
    end
    if addInterface
        fID=fopen(fileNameNew,'a');
        fprintf(fID,'%s\n',['*Element, type=COH2D4, elset=Interface']);
        dlmwrite(fileNameNew,CohEle_Matrix,'-append','precision',10);
        fclose(fID);
    end
    
    fID=fopen(fileNameNew,'a');
    for kpart=1:nPart
        fprintf(fID,'%s\n',['*Solid Section, elset=PART-',num2str(kpart),', material=PART']);
        fprintf(fID,'%s\n',[',']);
    end
    fprintf(fID,'%s\n',['*Solid Section, elset=MATRIX, material=MATRIX']);
    fprintf(fID,'%s\n',[',']);
    if YLmemb
        fprintf(fID,'%s\n',['*Solid Section, elset=YL, material=Surface']);
        fprintf(fID,'%s\n',['1.0']);
    end
    if addInterface
        fprintf(fID,'%s\n',['*Cohesive section, elset=Interface, material=matCOH']);
        fprintf(fID,'%s\n',['1.0']);
    end
    
    fprintf(fID,'%s\n',['*Nset, nset=NODEA']);
    fprintf(fID,'%s\n',[num2str(Node_A)]);
    fprintf(fID,'%s\n',['*Nset, nset=NODEB']);
    fprintf(fID,'%s\n',[num2str(Node_B)]);
    fprintf(fID,'%s\n',['*Nset, nset=NODEC']);
    fprintf(fID,'%s\n',[num2str(Node_C)]);
    fprintf(fID,'%s\n',['*Nset, nset=NODED']);
    fprintf(fID,'%s\n',[num2str(Node_D)]);
    fprintf(fID,'%s\n',['*Nset, nset=nOUT']);
    fprintf(fID,'%s\n',['nodeXX,nodeYY,nodeXY']);
    
    fprintf(fID,'%s\n',['*Nset, nset=NODEXN,unsorted']);
    fclose(fID);
    dlmwrite(fileNameNew,NodeXN,'-append','precision',10);
    fID=fopen(fileNameNew,'a');
    
    fprintf(fID,'%s\n',['*Nset, nset=NODEXP,unsorted']);
    fclose(fID);
    dlmwrite(fileNameNew,NodeXP,'-append','precision',10);
    fID=fopen(fileNameNew,'a');
    
    fprintf(fID,'%s\n',['*Nset, nset=NODEYN,unsorted']);
    fclose(fID);
    dlmwrite(fileNameNew,NodeYN,'-append','precision',10);
    fID=fopen(fileNameNew,'a');
    
    fprintf(fID,'%s\n',['*Nset, nset=NODEYP,unsorted']);
    fclose(fID);
    dlmwrite(fileNameNew,NodeYP,'-append','precision',10);
    fID=fopen(fileNameNew,'a');
    
    fprintf(fID,'%s \n',['*Equation']);
    fprintf(fID,'%s \n',['3']);
    fprintf(fID,'%s \n',['NODEB,1,1.0,NODEA,1,-1.0,nodeXX,1,-1.0']);
    fprintf(fID,'%s \n',['*Equation']);
    fprintf(fID,'%s \n',['2']);
    fprintf(fID,'%s \n',['NODEB,2,1.0,NODEA,2,-1.0']);
    
    fprintf(fID,'%s \n',['*Equation']);
    fprintf(fID,'%s \n',['4']);
    fprintf(fID,'%s \n',['NODEC,1,1.0,NODEA,1,-1.0,nodeXX,1,-1.0,nodeXY,1,-1.0']);
    fprintf(fID,'%s \n',['*Equation']);
    fprintf(fID,'%s \n',['3']);
    fprintf(fID,'%s \n',['NODEC,2,1.0,NODEA,2,-1.0,nodeYY,2,-1.0']);
    
    fprintf(fID,'%s \n',['*Equation']);
    fprintf(fID,'%s \n',['3']);
    fprintf(fID,'%s \n',['NODED,1,1.0,NODEA,1,-1.0,nodeXY,1,-1.0']);
    fprintf(fID,'%s \n',['*Equation']);
    fprintf(fID,'%s \n',['3']);
    fprintf(fID,'%s \n',['NODED,2,1.0,NODEA,2,-1.0,nodeYY,2,-1.0']);
    
    fprintf(fID,'%s \n',['*Equation']);
    fprintf(fID,'%s \n',['3']);
    fprintf(fID,'%s \n',['NODEXP,1,1.0,NODEXN,1,-1.0,nodeXX,1,-1.0']);
    fprintf(fID,'%s \n',['*Equation']);
    fprintf(fID,'%s \n',['2']);
    fprintf(fID,'%s \n',['NODEXP,2,1.0,NODEXN,2,-1.0']);
    
    fprintf(fID,'%s \n',['*Equation']);
    fprintf(fID,'%s \n',['3']);
    fprintf(fID,'%s \n',['NODEYP,1,1.0,NODEYN,1,-1.0,nodeXY,1,-1.0']);
    fprintf(fID,'%s \n',['*Equation']);
    fprintf(fID,'%s \n',['3']);
    fprintf(fID,'%s \n',['NODEYP,2,1.0,NODEYN,2,-1.0,nodeYY,2,-1.0']);
    
    fprintf(fID,'%s\n',['*End Part']);
    
     %
    % assembly
    %
    fprintf(fID,'%s\n',['**']);
    fprintf(fID,'%s\n',['**']);
    fprintf(fID,'%s\n',['** ASSEMBLY']);
    fprintf(fID,'%s\n',['**']);
    fprintf(fID,'%s\n',['*Assembly, name=Assembly']);
    fprintf(fID,'%s\n',['**']);
    fprintf(fID,'%s\n',['*Instance, name=Part-1-1, part=Part-1']);
    fprintf(fID,'%s\n',['*End Instance']);
    fprintf(fID,'%s\n',['**  ']);
    fprintf(fID,'%s\n',['*End Assembly']);
    fprintf(fID,'%s\n',['*Section controls, name=haqiCoh, viscosity=1.0e-3']);
     % Constrain the reference point DOF
    fprintf(fID,'%s \n',['*BOUNDARY']);
    fprintf(fID,'%s \n',['Part-1-1.nodeXX, 2, 2, 0']);
    fprintf(fID,'%s \n',['Part-1-1.nodeYY, 1, 1, 0']);
    fprintf(fID,'%s \n',['Part-1-1.nodeXY, 2, 2, 0']);
    fprintf(fID,'%s \n',['Part-1-1.',num2str(Node_Alpha),', 1, 2, 0']);
    fprintf(fID,'%s \n',['Part-1-1.',num2str(Node_Beta),', 2, 2, 0']);
    
    fprintf(fID,'%s\n',['**']);
    fprintf(fID,'%s\n',['** MATERIALS']);
    fprintf(fID,'%s\n',['**']);
    fprintf(fID,'%s\n',['*Material, name=PART']);
    fprintf(fID,'%s\n',['*Elastic']);
    fprintf(fID,'%s\n',['400000,0.26']);
    
    fprintf(fID,'%s\n',['*Material, name=MATRIX']);
    fprintf(fID,'%s\n',['*Elastic']);
    fprintf(fID,'%s\n',['70000,0.3']);
    if YLmemb
        fprintf(fID,'%s\n',['*Material, name=Surface']);
        fprintf(fID,'%s\n',['*Elastic']);
        fprintf(fID,'%s\n',['100000000000,0.1']);
    end
    if addInterface
        fprintf(fID,'%s\n','*Material, name=matCOH');
        fprintf(fID,'%s\n','*elastic, type=traction');
        fprintf(fID,'%s\n',['1.0E+8,1.0E+8,1.0E+8']);
        fprintf(fID,'%s\n','*damage initiation,CRITERION=QUADS');
        fprintf(fID,'%G %s %G %s %G\n',75, ',', 75, ',', 50);
        fprintf(fID,'%s\n','*Damage Evolution, type=DISPLACEMENT');
        fprintf(fID,'%s\n','0.088671,');
    end
    fprintf(fID,'%s\n',['** ----------------------------------------------------------------']);
        fprintf(fID,'%s\n',['** ']);
        fprintf(fID,'%s\n',['** STEP: Volume Calculation']);
        fprintf(fID,'%s\n',['** ']);
        fprintf(fID,'%s\n',['*Step,NAME=volumeCalculation, NLGEOM=NO, inc=1000 ']);
        fprintf(fID,'%s\n',['*STATIC ']);
        fprintf(fID,'%s\n',['1 , 1']);
        fprintf(fID,'%s\n',['*Output, field']);
        fprintf(fID,'%s\n',['*Element Output']);
        fprintf(fID,'%s\n',['EVOL']);
        fprintf(fID,'%s\n',['*End Step']);
        fprintf(fID,'%s\n',['** ']);
        fprintf(fID,'%s\n',['** STEP: Step-1']);
        fprintf(fID,'%s\n',['** ']);
        fprintf(fID,'%s\n',['*Step, name=Step-1, nlgeom=NO, convert sdi=NO,inc=100000']);
        fprintf(fID,'%s\n',['*Static']);
        fprintf(fID,'%s\n',['1 , 1']);

        fprintf(fID,'%s\n',['*CONTROLS,PARAMETERS=FIELD,FIELD=DISPLACEMENT']);
        fprintf(fID,'%s\n',['0.1,10.0,,,100.0,1.0,1.0,10.0 ']);
        fprintf(fID,'%s\n',['*controls,parameters=time incrementation ']);
        fprintf(fID,'%s\n',[' 8,10 , , , , , ,20']);
        fprintf(fID,'%s\n',['*CONTROLS,PARAMETERS=LINE SEARCH']);
        fprintf(fID,'%s\n',['30,,,,']);
        fprintf(fID,'%s \n',['*BOUNDARY']);
        fprintf(fID,'%s\n',['Part-1-1.nodeXX, 1, 1, ' num2str(0.01*2*xBound)]);
        
        fprintf(fID,'%s\n',['** OUTPUT REQUESTS']);
%         fprintf(fID,'%s\n',['*Output, field, time interval=5e-03']);
        fprintf(fID,'%s\n',['*Output, field']);
        fprintf(fID,'%s\n',['*Node Output']);
        fprintf(fID,'%s\n',['U']);
        fprintf(fID,'%s\n',['*Element Output']);
        fprintf(fID,'%s\n',['S']);
        fprintf(fID,'%s\n',['E']);
        
        if addInterface
            fprintf(fID,'%s\n',['*Element Output, elset=Part-1-1.INTERFACE']);
            fprintf(fID,'%s\n',['STATUS']);
            fprintf(fID,'%s\n',['SDEG']);
            fprintf(fID,'%s\n',['DMICRT']);
            fprintf(fID,'%s\n',['QUADSCRT']);
        end
        fprintf(fID,'%s\n',['*Output, history, frequency=1']);
        fprintf(fID,'%s\n',['*Node Output, nset=Part-1-1.nOut ']);
        fprintf(fID,'%s\n',['RF, U']);
        fprintf(fID,'%s\n',['*End Step']);
        disp('·ÖÎö¤·¤Þ¤¹£º')
        cd(dirname);
        system(['echo y | abaqus job=',fileNameNew((length(dirname)+2):end-4),' int']);
        cd ..;
end
toc;