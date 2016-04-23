coord=final_points(:,[1:2,4]);
dist=pdist(coord);
tree=linkage(dist, 'single'); %linkn based on shortest distance
clusters=cluster(tree,'maxclust',7);
