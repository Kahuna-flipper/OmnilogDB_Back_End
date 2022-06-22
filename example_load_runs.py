
from IPython import embed
import pandas as pd
import seaborn as sns
import pylab as plt
import numpy as np

import os
basedir = os.path.abspath(os.path.dirname(__file__))

from phenom.phenom_db import *

''' CONNECTING TO DATABASES '''
# MYSQL - for metadata
from config import SQLALCHEMY_DATABASE_URI, MONGO_CNX

from sqlalchemy import create_engine
engine = create_engine(SQLALCHEMY_DATABASE_URI)
from sqlalchemy.orm import sessionmaker
Session = sessionmaker(bind=engine)
session = Session()

# Mongo for growth data
import pymongo
# Making a Connection with MongoClient
myclient = pymongo.MongoClient(MONGO_CNX['client'])
# database
mydb = myclient[MONGO_CNX['db']]
# collection
mongo_table= mydb[MONGO_CNX['collection']]

''' END DB CONNECTIONS'''

''' location of exported biolog data '''
DATA_DIR = r'data/'

from phenom import io
from phenom import analysis


def main():

    load_sample_data('sample_biolog_data_metadata_curated.xlsx', overwrite=True)
    
    # show_db_stats()
    
    
def load_sample_data(filename, overwrite=False):
    '''
    Loads curated biolog data into the data basedir
    
    Set overwrite = True to delete all data in the database and re-create all tables
    
    To load additional data to an existing database set overwrite = False
    
    '''

    ''' 0. clean database and load plate info ONLY RUN IF YOU WANT TO RECREATE THE DATABASE '''
    if overwrite == True:
        io.clean_dbs(engine, mongo_table)
        io.load_plate_info(engine)
        
    ''' 1. Load in example data '''
    all_data = pd.read_excel(DATA_DIR + filename)
    
    
    ''' 2. parse data into runinfo and growth_data '''
    res = io.parse_data_for_db(all_data, engine)
    run_info = res['run_info']
    growth_data = res['growth_data']
    
    
    ''' 3 load run_info & studies'''
    io.add_runinfo(run_info, engine) # use this to add new data
    
    
    
    ''' 4. load growth_info'''
    io.load_growth_data(growth_data, mongo_table)

    ''' 5. calculate and load all aucs and fold changes (fc) '''
    
    from phenom import analysis 
    runs = session.query(Runinfo).all()
    for run in runs:
        plate = run.plate_id
        data = io.get_plate_data(plate, session, mongo_table)
        aucs = analysis.calc_aucs(data)
        fc = analysis.calc_fcs(data)
        growth_metrics = pd.merge(aucs[['auc','auc_ratio','well','plate_id']], fc[['fold_change','log2FC','well','plate_id']], left_on=['plate_id','well'], right_on=['plate_id','well'])
        io.load_aucs(growth_metrics, engine)

def show_db_stats():
    runinfo = session.query(Runinfo)
    runinfo = pd.read_sql(runinfo.statement, runinfo.session.bind)
    
    runinfo = runinfo[runinfo['setup_time']!='Setup Time']
    runinfo['setup_time'] = pd.to_datetime(runinfo['setup_time'])
    res = runinfo.groupby('setup_time').count()['plate_id']
    monthly = res.resample('M').sum()
    cumulative = monthly.cumsum()
    print(cumulative)
    cumulative.plot()
    plt.ylabel('Cumulative plates run')
    plt.title('SBRG Omnilog Plates used')
    
    plt.show()
    

if __name__=='__main__':
    main()