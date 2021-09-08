%Graph Transformation Matching
%Written By: Z. Rozita Ghassabi (z.r.ghassabi@gmail.com)

%=======================================
%Example test code
%=======================================
%points1=[10 10;11 50;25 25;32 8;43 40];
%points2=[10 14;12 48;25 28;32 7;48 38];
%GTM(points1,points2,2)
   
function [ Nodes1 Nodes2 ] = GTM(points1,points2,K)
    tic;
    distp1=ComputeDistanceMatrix(points1(:,1:2));
    distp2=ComputeDistanceMatrix(points2(:,1:2));
    
    %medianp1=mean(mean(distp1));
    %medianp2=mean(mean(distp2));
    medianp1=computeMedian(distp1);
    medianp2=computeMedian(distp2);
    
    %K=2;
    Q1=points1(:,1:2);
    Q2=points2(:,1:2);

    AP1=BuildMedianKNNGraph3(distp1,K,medianp1);
    AP2=BuildMedianKNNGraph3(distp2,K,medianp2);
    %AP1=tmpBuildMedianGraph(distp1,medianp1);
    %AP2=tmpBuildMedianGraph(distp2,medianp2);

    numIteration=0;
    while isequal(AP1,AP2)==0
        numIteration=numIteration+1;
        jOut=FindOutlier(AP1,AP2);
        %disp(['it:',num2str(numIteration),'->',num2str(jOut)]);
        if jOut~=0
            Q1(jOut,:)=[];
            Q2(jOut,:)=[];
            
            distp1(:,jOut)=[];
            distp1(jOut,:)=[];
            distp2(:,jOut)=[];
            distp2(jOut,:)=[];
            %distp1=ComputeDistanceMatrix(Q1);
            %distp2=ComputeDistanceMatrix(Q2);

            %medianp1=mean(mean(distp1));
            %medianp2=mean(mean(distp2));
            medianp1=computeMedian(distp1);
            medianp2=computeMedian(distp2);

            AP1=BuildMedianKNNGraph3(distp1,K,medianp1);
            AP2=BuildMedianKNNGraph3(distp2,K,medianp2);
            %AP1=tmpBuildMedianGraph(distp1,medianp1);
            %AP2=tmpBuildMedianGraph(distp2,medianp2);
        end
    end
    tmElapsed=toc;
    disp(['Iteration count:',num2str(numIteration)]);
    disp(['Time:',num2str(tmElapsed),' sec']);
    %=======================================
    %For return same points
    %=======================================
    conVertices=GetConnectedVertices(AP1);
    if conVertices~=0
        Q1=Q1(conVertices,:);
        Q2=Q2(conVertices,:);
        %sameNodes=[Q1 Q2];
         Nodes1=Q1;
         Nodes2=Q2;
    else
        sameNodes='[NO MATCH POINTS]';
    end
    
    %=======================================
    %For return isomorphism graph matrix
    %=======================================
    %Q1=RemoveDisconnectedVertices(AP1);
    %Q2=RemoveDisconnectedVertices(AP2);
    %sameNodes=[Q1];
    %=======================================
end

function [ med ] = computeMedian( distMatrix )
    nRows=size(distMatrix,1);
    temp=distMatrix(2,:);
    for i=2:nRows
        temp = [temp distMatrix(i,:)];
    end
    temp=sort(temp);
    n=size(temp,2);
    n=(n-mod(n,2))/2;
    med=temp(n+1);
end

function [ distMatrix ] = ComputeDistanceMatrix( points )
    pointsCount=size(points,1);
    distMatrix=zeros(pointsCount,pointsCount);
    
    for elemID=1:1:pointsCount 
        for elem2=1:1:pointsCount 
            distMatrix(elemID,elem2)=sqrt((points(elemID,1)-points(elem2,1))^2 + (points(elemID,2)-points(elem2,2))^2);
        end
    end
end

function [ KNNGraph ] = BuildMedianKNNGraph( distMatrix,K,medianValue )
    pointCount=size(distMatrix,1);
    KNNGraph=zeros(pointCount,pointCount);
    %temp1=zeros(pointCount);
    %temp2=zeros(pointCount);
    %maxDist=0;
    distMatrix(distMatrix>medianValue)=0;
    
    for i=1:1:pointCount 
        temp1=distMatrix(i,:);
        
        if temp1(:)==0
        else
            temp2=sort(temp1(temp1>0));
            if (size(temp2,2)>K)
                maxDist=temp2(K);   
                temp1(temp1>maxDist)=0;
                temp1(temp1>0)=1;
            elseif (size(temp2,2)<K)
                temp1(:)=0; 
            else
                temp1(temp1>0)=1;
            end
            
        end
        KNNGraph(i,:)=temp1;
    end
end

function [ KNNGraph ] = BuildMedianKNNGraph3( distMatrix,K,medianValue )
    pointCount=size(distMatrix,1);
    KNNGraph=zeros(pointCount,pointCount);
    bigNumber=10e06;
        
    tmpDistanceMatrix=distMatrix;
    tmpDistanceMatrix(tmpDistanceMatrix>medianValue | tmpDistanceMatrix==0 )=bigNumber;
    tmpDistanceMatrix=sort(tmpDistanceMatrix')';
    
    for i=1:1:pointCount
        if tmpDistanceMatrix(i,K)~=bigNumber
            temp1=tmpDistanceMatrix(i,1:K);
            for j=1:K
                item=temp1(j);
                itemIndex=find(distMatrix(i,:)==item);
                KNNGraph(i,itemIndex)=1;
                KNNGraph(itemIndex,i)=1;
            end
        end
    end
end

function [ KNNGraph ] = BuildMedianKNNGraph2( distMatrix,K,medianValue )
    pointCount=size(distMatrix,1);
    KNNGraph=zeros(pointCount,pointCount);
    %temp1=zeros(pointCount);
    %temp2=zeros(pointCount);
    %maxDist=0;
    distMatrix(distMatrix>medianValue)=0;
    tempCount(1:pointCount)=0;
    tempK=0;
    
    for i=1:1:pointCount 
        temp1=distMatrix(i,i:pointCount);
        tempK=K;%-tempCount(i);
        
        if tempK<=0
            continue;
        end
                    
        if temp1(:)==0
        else
            temp2=sort(temp1(temp1>0));
            if (size(temp2,2)>tempK)
                maxDist=temp2(tempK);   
                temp1(temp1>maxDist)=0;
                temp1(temp1>0)=1;
            elseif (size(temp2,2)<tempK)
                temp1(:)=0; 
            else
                temp1(temp1>0)=1;
            end
            
        end
        
        KNNGraph(i,i:pointCount)=temp1;
        KNNGraph(i:pointCount,i)=temp1';
        %tempCount=tempCount+KNNGraph(i,1:pointCount); 
        
        %for j=1:pointCount
        %    if KNNGraph(i,j)==1
        %        if j~=i
        %           tempCount(j)=tempCount(j)+1; 
        %        end
        %    end
        %end
    end
end

function [ KNNGraph ] = tmpBuildMedianGraph( distMatrix,medianValue )
    pointCount=size(distMatrix,1);
    KNNGraph=zeros(pointCount,pointCount);
    %temp1=zeros(pointCount);
    %temp2=zeros(pointCount);
    %maxDist=0;
    KNNGraph=distMatrix;
    KNNGraph(KNNGraph>medianValue)=0;
    KNNGraph(KNNGraph>0)=1;
end

function [ jOutlier ] = FindOutlier( matrixAP1,matrixAP2 )
    R=abs(matrixAP1-matrixAP2);
    tmp=sum(R);
    jOutlier=0;
    if sum(tmp)~=0
        nR=size(tmp,2);
        maxNum=max(tmp);
        for kMax=1:nR
            if tmp(kMax)==maxNum
                jOutlier=kMax;     
                break;
            end
        end
    end
    
end

function [ resNodesNumber ] = GetConnectedVertices( matrixAdj )
    tmp=sum(matrixAdj);
    n=0;
    sz=size(tmp,2);
    tmp1=0;
    for i=1:sz
    if tmp(i)~=0
        n=n+1;
        tmp1(n)=i;
    end    
    end
    resNodesNumber=tmp1;
end

function [ resMatrix ] = RemoveDisconnectedVertices( matrixAdj )
    tmp=sum(matrixAdj);
    n=0;
    sz=size(tmp,2);
    for i=1:sz
    if tmp(i)~=0
        n=n+1;
        tmp1(n)=i;
    end    
    end
    resMatrix=matrixAdj(tmp1,tmp1);
end

