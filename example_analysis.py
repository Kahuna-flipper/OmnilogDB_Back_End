from IPython import embed
import pandas as pd
import seaborn as sns
import pylab as plt
import numpy as np

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


from phenom import io
from phenom import analysis
from phenom import reports

def main():
    study_name = 'P. putida JBEI'
    matrix = analysis.get_study_auc_matrix(session, study_name, index=['strain','media'], param='auc_ratio')
    aucs = analysis.get_study_processed_data(session, study_name)
    
    # NOTE this analysis will create a PDF with figures of growth curves and bar plots for comparisons, it can take ~20 minutes to run depending on the number of samples being analyzed
    # Be sure that a folder named "reports" is in the base repo directory, its contents are ignored by git
    
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt
    
    # Create the PdfPages object to which we will save the pages:
    # The with statement makes sure that the PdfPages object is closed properly at
    # the end of the block, even if an Exception occurs.
    aucs['strain_media'] = aucs['strain']+' '+aucs['media']
    study_name_file = study_name.replace('/','-')
    with PdfPages('reports/' + study_name_file + '.pdf') as pdf:
        for plate_type in aucs['plate_type'].unique():
            print(plate_type)
            aucs2 = aucs[aucs['plate_type']==plate_type]
            print(len(aucs2))
            reports.PCA_plots(aucs2, pdf, value='auc', name=plate_type + ' AUC', columns=['strain_media'])
            
            reports.barcharts(aucs2, pdf, value='auc', columns=['strain','media'])
            
        
        data = analysis.get_study_raw_data(session, study_name, mongo_table)
        data['strain_media'] = data['strain']+' '+data['media']
        
        reports.ts_plots(data, pdf, group='strain_media')
    
    embed()
    
    




if __name__=='__main__':
    main()