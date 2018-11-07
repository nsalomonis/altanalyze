### SQL Interface

import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies

import time
import random
import math
import sqlite3
import export

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def filepath(filename):
    try:
        import unique ### local to AltAnalyze
        fn = unique.filepath(filename)
    except Exception:
        ### Should work fine when run as a script with this (AltAnalyze code is specific for packaging with AltAnalyze)
        dir=os.path.dirname(dirfile.__file__)
        try: dir_list = os.listdir(filename); fn = filename ### test to see if the path can be found (then it is the full path)
        except Exception: fn=os.path.join(dir,filename)
    return fn
    
##### SQLite Database Access ######

def createSchemaTextFile(species,platform,schema_text,DBname):
    schema_filename = filepath('AltDatabase/'+species+'/'+platform+'/'+DBname+'_schema.sql')
    export_data = export.ExportFile(schema_filename)

    ### We will need to augment the database with protein feature annotations for 
    export_data.write(schema_text)
    export_data.close()
   
def populateSQLite(species,platform,DBname,schema_text=None):
    global conn
    """ Since we wish to work with only one gene at a time which can be associated with a lot of data
    it would be more memory efficient transfer this data to a propper relational database for each query """
    
    db_filename = filepath('AltDatabase/'+species+'/'+platform+'/'+DBname+'.db') ### store in user directory
    schema_filename = filepath('AltDatabase/'+species+'/'+platform+'/'+DBname+'_schema.sql')
    
    ### Check to see if the database exists already and if not creat it
    db_is_new = not os.path.exists(db_filename)

    with sqlite3.connect(db_filename) as conn:
        if db_is_new:
            createSchemaTextFile(species,platform,schema_text,DBname)
            print 'Creating schema'
            with open(schema_filename, 'rt') as f:
                schema = f.read()
            #print schema
            conn.executescript(schema)
            
        else:
            print 'Database exists, assume schema does too.'
            #sys.exit()
            
    return conn ### User must now add data to the empty SQLite database
    
def connectToDB(species,platform,DBname):
    db_filename = filepath('AltDatabase/'+species+'/'+platform+'/'+DBname+'.db') ### store in user directory
    with sqlite3.connect(db_filename) as conn:
        return conn

def retreiveDatabaseFields(conn,ids,query):
    """ Retreive data from specific fields from the database """
    cursor = conn.cursor()
    
    #id = 'ENSG00000114127'
    #query = "select id, name, description, chr, strand from genes where id = ?"
    
    cursor.execute(query,ids) ### In this way, don't have to use %s and specify type
    
    ls=[]
    for row in cursor.fetchall():
        #id, name, description, chr, strand = row
        #print '%s %s %s %s %s' % (id, name, description, chr, strand)
        ls.append(row)
    return ls
        
def bulkLoading():
    import csv
    import sqlite3
    import sys
    
    db_filename = 'todo.db'
    data_filename = sys.argv[1]
    
    SQL = """insert into task (details, priority, status, deadline, project)
             values (:details, :priority, 'active', :deadline, :project)
          """
    
    with open(data_filename, 'rt') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        
        with sqlite3.connect(db_filename) as conn:
            cursor = conn.cursor()
            cursor.executemany(SQL, csv_reader)
        