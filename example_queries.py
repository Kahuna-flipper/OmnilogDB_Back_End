from IPython import embed
import pandas as pd
import seaborn as sns
import pylab as plt
import numpy as np

from phenom.phenom_db import *

import os
basedir = os.path.abspath(os.path.dirname(__file__))


''' CONNECTING TO DATABASES '''
from config import SQLALCHEMY_DATABASE_URI, MONGO_CNX
from sqlalchemy import create_engine
engine = create_engine(SQLALCHEMY_DATABASE_URI)
from sqlalchemy.orm import sessionmaker
Session = sessionmaker(bind=engine)
session = Session()


''' Query all runs '''
runs = session.query(Runinfo)
runinfo = pd.read_sql(runs.statement, runs.session.bind)
print(runinfo)


''' GETTING STRAIN GROWTH INFO '''
strain = '700447'  # an example strain to query for
runs = session.query(Runinfo).filter(Runinfo.strain==strain)
runinfo = pd.read_sql(runs.statement, runs.session.bind)

aucs = session.query(Growth).filter(Growth.plate_id.in_(runinfo['plate_id'].tolist()))
aucs = pd.read_sql(aucs.statement, aucs.session.bind)
print(aucs)


''' GETTING COMPOUND GROWTH INFO '''
compound = "2% NaCl"
metric = 'auc'
info = session.query(Plateinfo).filter(Plateinfo.name.contains(compound))
info = pd.read_sql(info.statement, info.session.bind)
    
aucs = session.query(Growth).filter(Growth.well.in_(info['well'].tolist()))
aucs = pd.read_sql(aucs.statement, aucs.session.bind)
runs = session.query(Runinfo)
runinfo = pd.read_sql(runs.statement, runs.session.bind)
aucs = pd.merge(aucs, runinfo)
grouped = aucs.groupby(['species','strain','media']).mean()[metric]
print(grouped)
