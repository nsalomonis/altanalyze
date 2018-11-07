import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies

from sklearn.metrics.cluster import adjusted_rand_score
import numpy as np

def ari(truelabel,predlabel):
    lab={}
    truelab=[]
    predlab=[]
    for line in open(truelabel,'rU').xreadlines():
        data = line.rstrip()
        t = string.split(data,'\t')
        lab[t[0]]=[int(t[1]),]
    for line in open(predlabel,'rU').xreadlines():
        data = line.rstrip()
        t = string.split(data,'\t')
        try:lab[t[0]].append(int(t[1]))
        except Exception: print "Sample missing true label"
    for key in lab:

        try:
            predlab.append(lab[key][1])
            truelab.append(lab[key][0])
        except Exception:
            print "Sample missing predicted label"
            continue
        
    print len(truelab)
    truelab=np.array(truelab)
    predlab=np.array(predlab)
  
    ari=adjusted_rand_score(truelab,predlab)
    return ari

#truelabel="/Volumes/Pass/Archive_Zeisel/SVMOutputs/groups.round1SVC_Results_max.txt"
#predlabel="/Volumes/Pass/Singlecellbest/Zeisel_upd/SVMOutputs/round1SVC_Results.txt"
#predlabel="/Volumes/Pass/Singlecellbest/Zeisel_upd/SVMOutputs/round1SVC_Results.txt"
#truelabel="/Volumes/Pass/Singlecellbest/Pollen_upd/SVMOutputs/groups.round1SVC_Results_max.txt"
#predlabel="/Volumes/Pass/Singlecellbest/Pollen_upd/SVMOutputs/round1SVC_Results.txt"
#predlabel="/Volumes/Pass/Data/Pollen_cluster.txt"
#predlabel="/Users/meenakshi/Usoskin_Sc3_test.txt"
#truelabel="/Volumes/Pass/Singlecellbest/Usoskin_upd/SVMOutputs/groups.round1SVC_Results_max.txt"
#predlabel="/Users/meenakshi/Downloads/k-11-Usoskin.txt"
#predlabel="/Users/meenakshi/Documents/ZeiselCluster.txt"
#truelabel="/Users/meenakshi/Desktop/groups.Pollen.txt"
#predlabel="/Users/meenakshi/Downloads/SC3_pollen.txt"
#predlabel="/Users/meenakshi/groups-filtered.txt"
truelabel="/Users/meenakshi/Documents/Singlecelldata/groups.CD34.v5_test.txt"
#predlabel="/Volumes/Pass/HCA_ICGS/ICGS_complete/SVMOutputs/round1SVC_Results_old30.txt"
#truelabel="/Volumes/Pass/HCA_ICGS/ICGS_complete/SVMOutputs/round1SVC_Results.txt"
#truelabel="/Volumes/Pass/HCA_Umap/SVMOutputs/groups.round1SVC_Results_max3.txt"
#truelabel="/Users/meenakshi/cluster-filtered.txt"
#predlabel="/Volumes/Pass/Immune-complete/SVMOutputs/round1SVC_Results_HGV_UMAP_PR_0.2_0.2.txt"
#truelabel="/Volumes/Pass/paperResults/tsne_in/Known_groups.txt"
predlabel="/Volumes/Pass/ICGS2_testrun/ICGS-NMF/FinalGrous_Merged.txt"
#truelabel="/Volumes/Pass/Final_scicgs/groups.round1SVC_Results-orig.txt"
#predlabel="/Volumes/Pass/HCA_Umap/SVMOutputs/round1SVC_Results.txt"
#predlabel="/Volumes/Pass/Deng/SVMOutputs/round1SVC_Results.txt"
#truelabel="/Volumes/Pass/Deng/groups.Deng.txt"
#truelabel="/Users/meenakshi/Documents/paperResults/tsne_in/donor_all.txt"
#predlabel="/Users/meenakshi/Documents/paperResults/tsne_in/round1SVC_Results_HGV_UMAP_PR_0.2_0.2.txt"
arival=ari(truelabel,predlabel)
print str(arival)