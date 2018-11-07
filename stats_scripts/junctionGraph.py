import os, sys, string

def createFeatures(junction_coordinates):
    features={}
    for junction_id in junction_coordinates:
        location = junction_coordinates[junction_id]
        pos1,pos2 = string.split(string.split(location,':')[-1],'-')
        pos = [int(pos1),int(pos2)]
        pos.sort()
        features[junction_id] = pos
    return features

def createFeaturesFromEvents(gene_db):
    event_cluster_db={}
    cluster_name='NewClu_1'
    cluster_id=1
    count=1
    for gene in gene_db:
        features={}
        junction_to_event={}
        event_coordinates = gene_db[gene]
        for event in event_coordinates:
            
            junction1, junction2 = string.split(event,'|')
            junction1 = string.join(string.split(junction1,':')[1:],':')
            coordinates1,coordinates2 = string.split(event_coordinates[event],'|')
                        
            pos1,pos2 = string.split(string.split(coordinates1,':')[-1],'-')
            pos = [int(pos1),int(pos2)]
            pos.sort(); features[junction1] = pos
        
            pos1,pos2 = string.split(string.split(coordinates2,':')[-1],'-')
            pos = [int(pos1),int(pos2)]
            pos.sort(); features[junction2] = pos
            try: junction_to_event[junction1].append(event)
            except Exception: junction_to_event[junction1]=[event]
            try: junction_to_event[junction2].append(event)
            except Exception: junction_to_event[junction2]=[event]
        cluster_junction,cluster_name,cluster_id,count = filterByLocalJunctionExp(features,cluster_name,cluster_id,count)

        for junction in cluster_junction:
            events = junction_to_event[junction]
            for event in events:
                event_cluster_db[event]=cluster_junction[junction]
    return event_cluster_db

def filterByLocalJunctionExp(features,cluster_name,cluster_id,count):
    junctions_to_compare={}
    overlapping_junctions_exp={}
    ovelapping_pos={}
    existing=[]
    overlapping_junctions_test={}
    for feature in features:
            pos1,pos2 = features[feature]
            for f2 in features:
                    flag=False
                    if f2!=feature:
                        alt_pos1,alt_pos2 = features[f2]
                        positions = [pos1,pos2,alt_pos1,alt_pos2]
                        positions.sort()
                        diff = positions.index(pos2)-positions.index(pos1)
                        if diff!=1: ### Hence the two junctions are overlapping
                            flag=True
                        else:
                            diff = positions.index(alt_pos2)-positions.index(alt_pos1)
                            if diff!=1:
                                flag=True ### Hence the two junctions are overlapping
                                    
                        if flag==True:
                            if feature not in existing and f2 not in existing:
                                count=count+1
                                overlapping_junctions_test[count]=[feature,]
                                overlapping_junctions_test[count].append(f2)
                                existing.append(feature)
                                existing.append(f2)
                                    
                            elif feature in existing and f2 not in existing:
                                for i in overlapping_junctions_test:
                                    if feature in overlapping_junctions_test[i]:
                                        overlapping_junctions_test[i].append(f2)
                                        existing.append(f2)
                            elif f2 in existing and feature not in existing:
                                for i in overlapping_junctions_test:
                                    if f2 in overlapping_junctions_test[i]:
                                        overlapping_junctions_test[i].append(feature)
                                        existing.append(feature)   
                            elif feature in existing and f2 in existing:
                                for i in overlapping_junctions_test:
                                    if feature in overlapping_junctions_test[i]:
                                            loc1=i
                                    if f2 in overlapping_junctions_test[i]:
                                        loc2=i
                                if loc1!=loc2:
                                    for jun in overlapping_junctions_test[loc2]:
                                        if jun not in overlapping_junctions_test[loc1]:
                                            overlapping_junctions_test[loc1].append(jun)
                                    del overlapping_junctions_test[loc2]

    cluster_junction={}
    #Finding clusters and corresponding junctions
    for count in overlapping_junctions_test:
        for feature in overlapping_junctions_test[count]:
                cluster_junction[feature]=cluster_name
                #print feature,cluster_name
                pass
        cluster_id+=1
        cluster_name=string.join('NewClu_'+str(cluster_id))
        cluster_name=cluster_name.replace(" ","")
    return cluster_junction,cluster_name,cluster_id,count
            
if __name__ == '__main__':
        
    uids = {1:'chr8:134583936-134558169',
    2:'chr8:134583936-134574921',
    3:'chr8:134558017-134519564',
    4:'chr8:134558017-134511432',
    5:'chr8:134511378-134478333',
    6:'chr8:134478137-134477336',
    7:'chr8:134475657-134472180',
    8:'chr8:134474118-134472180',
    9:'chr8:134583936-134579563',
    10:'chr8:134511378-134488843',
    11:'chr8:134487962-134478333',
    12:'chr8:134488521-134488317',
    13:'chr8:134478137-134477200',
    14:'chr8:134475657-134474237'}

    cluster_name='NewClu_1'
    cluster_id=1
    count=1
        
    features = createFeatures(uids)
    cluster_junction,cluster_name,cluster_id,count = filterByLocalJunctionExp(features,cluster_name,cluster_id,count)
    cluster_junction,cluster_name,cluster_id,count = filterByLocalJunctionExp(features,cluster_name,cluster_id,count)