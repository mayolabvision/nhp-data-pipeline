function [x,y,c1,c2] = rsc_bin_vectors(Rtbl,Tclust,binID)

% attach binIDs to cluster table
Tclust.binID = binID;

% build lookup keys
key_clust = strcat(string(Tclust.sess_name),"_",string(Tclust.probe_index),"_",string(Tclust.cluster_id));
key_n1 = strcat(string(Rtbl.sess_name),"_",string(Rtbl.probe_index),"_",string(Rtbl.n1_clust));
key_n2 = strcat(string(Rtbl.sess_name),"_",string(Rtbl.probe_index),"_",string(Rtbl.n2_clust));

% find matching clusters
[~,loc1] = ismember(key_n1,key_clust);
[~,loc2] = ismember(key_n2,key_clust);

% bin IDs for each neuron in the pair
x = Tclust.binID(loc1);
y = Tclust.binID(loc2);

% rsc values
c1 = vertcat(Rtbl.rsc{:});

% rsc values
c2 = vertcat(Rtbl.rsig{:});

end