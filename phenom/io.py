import pandas as pd
from .phenom_db import *

ABX_PLATES = ['PM11C','PM12B']  
MY_PLATES = ['PM01','PM02','PM02A','PM03B','PM04A','PM09','PM10','PM11C','PM12B'] # only search for plates that we have run
METADATA = ['study', 'Plate Type','species','strain','media','rep', 'plate_id','Setup Time']
DB_METADATA = ['study', 'plate_type','species','strain','media','rep', 'plate_id','setup_time']
WELLS = [u'A01', u'A02', u'A03', u'A04', u'A05', u'A06', u'A07', u'A08', u'A09',
       u'A10', u'A11', u'A12', u'B01', u'B02', u'B03', u'B04', u'B05', u'B06',
       u'B07', u'B08', u'B09', u'B10', u'B11', u'B12', u'C01', u'C02', u'C03',
       u'C04', u'C05', u'C06', u'C07', u'C08', u'C09', u'C10', u'C11', u'C12',
       u'D01', u'D02', u'D03', u'D04', u'D05', u'D06', u'D07', u'D08', u'D09',
       u'D10', u'D11', u'D12', u'E01', u'E02', u'E03', u'E04', u'E05', u'E06',
       u'E07', u'E08', u'E09', u'E10', u'E11', u'E12', u'F01', u'F02', u'F03',
       u'F04', u'F05', u'F06', u'F07', u'F08', u'F09', u'F10', u'F11', u'F12',
       u'G01', u'G02', u'G03', u'G04', u'G05', u'G06', u'G07', u'G08', u'G09',
       u'G10', u'G11', u'G12', u'H01', u'H02', u'H03', u'H04', u'H05', u'H06',
       u'H07', u'H08', u'H09', u'H10', u'H11', u'H12']
       
       
def parse_data_for_db(data, engine):
    '''
    Takes in raw data from biolog exports and parses it into runinfo and growth data tables
    '''
    from sqlalchemy.orm import sessionmaker
    Session = sessionmaker(bind=engine)
    session = Session()
    
    # generating unique plate IDs - these are unique by date of run and position in the machine
    tmp=data[['Setup Time','Position']].drop_duplicates() 
    tmp2=data[['Setup Time','Position','study']].drop_duplicates() 
    print(tmp2.groupby('study').count())

    last_plate = session.query(Runinfo).order_by(Runinfo.plate_id.desc()).first() # get the latest plate
    if last_plate==None: # if the dB is empty:
        tmp['plate_id']=range(1,len(tmp)+1) # for interger plate IDs
    else:
        next_plate = last_plate.plate_id+1
        tmp['plate_id']=range(next_plate,len(tmp)+next_plate)
    
    data = pd.merge(data, tmp, left_on=['Setup Time','Position'], right_on=['Setup Time','Position']) 
    
    data = data[METADATA+WELLS+['Hour']]
    
    data = data.rename(columns={
        'Plate Type':'plate_type',
        'Hour':'hour',
        'Setup Time':'setup_time',
        'Position':'position'
    })
    
    runinfo = data[DB_METADATA].drop_duplicates()
    growth_data = data[['plate_id','hour'] + WELLS]
    
    return {'run_info':runinfo, 'growth_data':growth_data}
    
    
def add_runinfo(runinfo, engine):
    
    studies = pd.DataFrame(runinfo['study'].drop_duplicates())
    
    from sqlalchemy.orm import sessionmaker
    Session = sessionmaker(bind=engine)
    session = Session()
    last_study = session.query(Studies).order_by(Studies.study_id.desc()).first() # get the latest study
    last_plate = session.query(Runinfo).order_by(Runinfo.plate_id.desc()).first() # get the latest plate
    
    if last_study==None:
        next_study = 0
    else:
        next_study = last_study.study_id+1
    
    studies['study_id'] = range(next_study,next_study+len(studies))
        
    studies2runinfo = runinfo[['study','plate_id']].drop_duplicates()
    studies2runinfo = pd.merge(studies, studies2runinfo, left_on='study', right_on='study')
    del studies2runinfo['study']
    studies = studies.rename(columns={'study':'study_name'})
    
    runinfo.to_sql('runinfo',engine, if_exists='append',index=False)
    studies.to_sql('studies',engine, if_exists='append',index=False)
    studies2runinfo.to_sql('study2runinfo',engine, if_exists='append',index=False)
    
def load_aucs(aucs, engine):
    for plate in aucs['plate_id'].unique():
        print('loading AUCs for plate %i'%plate)
        plate_aucs = aucs[aucs['plate_id']==plate]
        plate_aucs.to_sql('growth',engine, if_exists='append',index=False)
    
    
def load_growth_data(growth_data, mongo_table):
    for plate in growth_data['plate_id'].unique():
        to_load = growth_data[growth_data['plate_id']==plate]
        print('loading plate_id: %i'%plate)
        to_load = to_load.astype({'hour':float})
        to_load.index = to_load['hour']
        plate_data_dict = to_load.to_dict("records")
        mongo_table.insert_one({"index":str(plate),"data":plate_data_dict})
    

def clean_dbs(engine, mongo_table):
    from sqlalchemy.ext.declarative import declarative_base
    from .phenom_db import Growth, Pca, Plateinfo, Study2runinfo, Studies, Runinfo, Base
    Base.metadata.drop_all(engine)
    Base.metadata.create_all(engine)
    
    mongo_table.drop()
    
def load_plate_info(engine):
    # loads all plate info (e.g. well compounds, cas #, etc to plateinfo table
    from glob import glob
    plate_files = glob('data\\plate_layouts\\*.csv')
    print('loading plate info for %s'%(', '.join(MY_PLATES)))
    all_plate_info = pd.DataFrame()
    for plate_file in plate_files:
        plate = plate_file.split('\\')[-1].replace('.csv','')
        if plate in MY_PLATES:
            info = parse_plate_info(plate)
            all_plate_info = pd.concat([all_plate_info, info])
            
    new_rows = len(all_plate_info)
    plates_to_add = len(all_plate_info['plate'].unique())
    if new_rows/plates_to_add!=96:
        print('are your sure all plate info was added?')
        embed()
    print('\n\n================================')
    print('Adding %i info rows for %i plates to the database'%(new_rows, plates_to_add))
    print('================================\n\n')
    
    all_plate_info.to_sql('plateinfo',engine, if_exists='append',index=False)
            
            
def parse_plate_info(plate):

    columns = pd.read_csv('data\\plate_layouts\\%s.csv'%plate, header=None, encoding='latin1') # encoding must be latin1 to overcome issues with CAS column
    columns = columns.rename(columns={0:'plate',1:'well',2:'name',3:'type', 8:'kegg',9:'cas'})
    columns = columns[['plate','well','name','type','kegg','cas']]
    
    
    if plate in ABX_PLATES: # ABX plates summarize data - need to break it out
        abx = ['1x','2x','3x','4x']
        out = []
        count=0
        for c in WELLS:
            name = columns.loc[columns.index[int(count/4)],'name']+'_'+abx[count%4]
            type = columns.loc[columns.index[int(count/4)],'type']
            kegg = columns.loc[columns.index[int(count/4)],'kegg']
            cas = columns.loc[columns.index[int(count/4)],'cas']
            out.append({'plate':plate, 'well':c, 'name':name,'type':type,'kegg':kegg,'cas':cas})
            count+=1
        columns = pd.DataFrame(out)
        
    return columns
    
    
def get_all_mongo_data(session, mongo_table):
    '''
    Exports all data from sql and nosql to a folder with files:
    runinfo.xlsx = metadata for all runs
    all data in csvs by plateid
    
    '''
    
    metadata = session.query(Runinfo)
    metadata = pd.read_sql(metadata.statement, metadata.session.bind)
    
    runinfo = session.query(Runinfo)
    runinfo = pd.read_sql(runinfo.statement, runinfo.session.bind)
    print('collecting growth data for %i plates'%len(runinfo))
    out = pd.DataFrame()
    for i in runinfo.index:
        plate = runinfo.loc[i,'plate_id']
        info = runinfo.loc[i,:]
        data_from_db = mongo_table.find_one({"index":str(plate)})
        df = pd.DataFrame(data_from_db["data"])
        out = pd.concat([out,df])
    
    out = pd.merge(metadata, out, left_on='plate_id', right_on='plate_id')
    
    
    return out
    
def get_growth_data(runinfo, session, mongo_table):
    out = pd.DataFrame()
    for i in runinfo.index:
        plate = runinfo.loc[i,'plate_id']
        info = runinfo.loc[i,:]
        data_from_db = mongo_table.find_one({"index":str(plate)})
        df = pd.DataFrame(data_from_db["data"])
        out = pd.concat([out,df])
    
    out = pd.merge(runinfo, out, left_on='plate_id', right_on='plate_id')
    
    return out
    
def get_plate_data(plate_id, session, mongo_table):
    runinfo = session.query(Runinfo).filter(Runinfo.plate_id==plate_id)
    runinfo = pd.read_sql(runinfo.statement, runinfo.session.bind)
    data_from_db = mongo_table.find_one({"index":str(plate_id)})
    df = pd.DataFrame(data_from_db["data"])
    
    out = pd.merge(runinfo, df, left_on='plate_id', right_on='plate_id')
    return out
    
    
def export_db_data_to_folder():
    '''
    Exports all data from sql and nosql to a folder with files:
    runinfo.xlsx = metadata for all runs
    all data in csvs by plateid
    
    '''
    
    runinfo = session.query(Runinfo)
    runinfo = pd.read_sql(runinfo.statement, runinfo.session.bind)
    
    runinfo.to_excel('exports/runinfo.xlsx')
    
    for i in runinfo.index:
        plate = runinfo.loc[i,'plate_id']
        info = runinfo.loc[:,i]
        data_from_db = mongo_table.find_one({"index":str(plate)})
        df = pd.DataFrame(data_from_db["data"])
        print('exporting plate %i'%plate)
        writer = pd.ExcelWriter('exports/plate_%i.xlsx'%plate)
        info.to_excel(writer,'metadata')
        df.to_excel(writer,'data')
        writer.save()
        writer.close()
    