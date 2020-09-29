function plotASCAcontrast(results,alpha)

% Plots the hierachical clusterring of the mean levels within the factor of
% interest, and puts on the inferential statistics onto the graph. 
% at level alpha (default = 0.05) the line is emphaised. 

if nargin==1
    alpha = 0.05;
end
% draw dendrogram and put in labels
dendrogram(results.linkageZ);
hold on
c = 0;
for i=1:size(results.linkageZ,1)-1
    yax = mean(results.linkageZ(i:i+1,3));
    if results.hierachical_contrasts(i,3)<alpha
        c = c+1;
    end
    if c==1
        yline(yax,'color','red', 'linewidth',3);
    end
    yline(yax,'color','red', 'linewidth',1);
    
    text(0.1,yax,['p = ' num2str(results.hierachical_contrasts(i,3))])
end
hold off;

