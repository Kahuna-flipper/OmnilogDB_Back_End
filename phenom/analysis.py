
from IPython import embed
import pandas as pd
import seaborn as sns
import pylab as plt
import numpy as np
from .phenom_db import *
from scipy.optimize import curve_fit
from sklearn.cluster import KMeans
from scipy.signal import savgol_filter
import scipy.stats

ABX_PLATES = ['PM11C','PM12B']  
METADATA = ['study', 'plate_type','species','strain','media','rep', 'hour','plate_id']
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


            
def calc_plate_fc(data, plate, interval=5):
    ''' 
    takes plate data in well (cols) x hour (rows) format and calculates fold change to last timepoint (or range of fimepoints
    specified by interval
    returns all fcs, log2FC as well as fcs relative neg
    '''
    
    columns = pd.read_csv('data/plate_layouts/%s.csv'%plate, header=None)
    columns = columns.rename(columns={0:'plate',1:'well',2:'name',3:'type', 8:'kegg',9:'cas'})
    
    out = []
    count=0
    abx = ['1x','2x','3x','4x']
    for c in WELLS:
        data2 = data[['Hour', c]]
        data2 = data2.astype({c:float})
        
        if plate in ABX_PLATES:
            name = columns.loc[columns.index[count/4],'name']+'_'+abx[count%4]
            type = columns.loc[columns.index[count/4],'type']
            kegg = columns.loc[columns.index[count/4],'kegg']
            cas = columns.loc[columns.index[count/4],'cas']
        else:
            idx = columns[columns['well']==c].index[0]
            name = columns.loc[idx,'name']
            type = columns.loc[idx,'type']
            kegg = columns.loc[idx,'kegg']
            cas = columns.loc[idx,'cas']
        
        
        out.append({'well':c, 'name':name,'type':type,'kegg':kegg,'cas':cas})
        count+=1
    columns = pd.DataFrame(out)
    
    data = data.sort_values('Hour')
    for c in WELLS:
        data = data.astype({c:float})
    
    #interval fold_change (from end time points to start time points)
    start = data[:interval].mean()
    end = data[-interval:].mean()
    fc = end/start
    fc.name='fold_change'
    fc = pd.DataFrame(fc)
    fc = fc.loc[WELLS,:]
    fc['log2FC'] = np.log2(fc)
    
    out = pd.merge(fc, columns[['well','name','type','kegg','cas']], left_index=True, right_on='well')
    return out
    
    
def calc_plate_kinetics(data, plate):

    out = pd.DataFrame(index=WELLS)
    functions = ['GompertzFunction','RichardsFunction','LogisticFunction']
    max_value = []
    area = []
    avg_height = []
    lag_time = []
    slope = []

    for well in WELLS:
        fit = np.zeros(3)
        growth = np.array(data[well])
        time = np.array(data['Hour'])
        growth = savgol_filter(growth, 11, 2) # window size 51, polynomial order 3

        while True:
            try:
                parameters1,pop = curve_fit(GompertzFunction,time,growth,maxfev=5000)
            except:
                #growth = savgol_filter(growth, 51, 2) # window size 51, polynomial order 2
                parameters1 = [1,1,1]
            break

        while True:
            try:
                parameters2,pop = curve_fit(RichardsFunction,time,growth,maxfev=5000)
            except:
                #growth = savgol_filter(growth, 51, 2) # window size 51, polynomial order 2
                parameters2=[0.1,0.1,0.1,0.1,0.1]
            break


        while True:
            try:
                parameters3,pop = curve_fit(LogisticFunction,time,growth,maxfev=5000)
            except:
                #growth = savgol_filter(growth, 51, 2) # window size 51, polynomial order 2
                parameters3=[1,1,1]
            break
        

        fit[0] = np.sum((GompertzFunction(time,*parameters1)-growth)**2)
        fit[1] = np.sum((RichardsFunction(time,*parameters2)-growth)**2)
        fit[2] = np.sum((LogisticFunction(time,*parameters3)-growth)**2)
        best_fit = np.argmin(fit)
    
        if(best_fit==0):
            fit_data = GompertzFunction(time,*parameters1)
            
        if(best_fit==1):
            fit_data = RichardsFunction(time,*parameters2)
            
        if(best_fit==2):
            fit_data = LogisticFunction(time,*parameters3)

        temp_max_value = np.max(fit_data)
        temp_area = np.trapz(fit_data)
        temp_avg_height = np.mean(fit_data)
        derivative_signal = np.zeros([np.size(fit_data),2])
        derivative_signal[:,0] = time

        for i in range(1,np.shape(derivative_signal)[0]):
            derivative_signal[i,1] = (fit_data[i]-fit_data[i-1])/(time[i]-time[i-1])
            
        temp_lag_time = time[np.argmax(derivative_signal[:,1])]
        temp_slope = np.max(derivative_signal)

        max_value.append(temp_max_value)
        area.append(temp_area)
        avg_height.append(temp_avg_height)
        lag_time.append(temp_lag_time)
        slope.append(temp_slope)

    out['Max Value'] = max_value
    out['AUC'] = area
    out['Avg Height'] = avg_height
    out['Lag time'] = lag_time
    out['Slope'] = slope
    
    return out

def calc_activity_index(param_set,plate,k):

    #param_set  = calc_plate_kinetics(data, plate)
    kmeans = KMeans(n_clusters=k, random_state=0).fit(param_set)
    param_set['Labels'] = kmeans.labels_

    avg_areas = []

    for label in np.linspace(0,k-1,k):
        temp_data = param_set[param_set['Labels']==label]
        avg_areas.append(np.average(temp_data['AUC']))
    
    activity_index_dataframe = pd.DataFrame(avg_areas,index=np.linspace(0,k-1,k),columns=['Avg AUC'])
    activity_index_dataframe = activity_index_dataframe.sort_values(by=['Avg AUC'],ascending=True)
    AI = np.linspace(0,k-1,k)
    activity_index_dataframe['Activity Index'] = AI

    ai_for_data = []

    for well in WELLS:
        label = param_set.loc[well]['Labels']
        AI = activity_index_dataframe.loc[label]['Activity Index']
        ai_for_data.append(AI)

    param_set['Activity Index'] = ai_for_data

    return param_set

            
            
def calc_plate_auc(data, plate):
    ''' 
    takes plate data in well (cols) x hour (rows) format and calculates auc for each curve
    returns all aucs as well as auc/auc of neg
    '''
    out = []
    
    columns = pd.read_csv('data/plate_layouts/%s.csv'%plate, header=None)
    columns = columns.rename(columns={0:'plate',1:'well',2:'name',3:'type'})
    
    count=0
    abx = ['1x','2x','3x','4x']
    for c in WELLS:
        data2 = data[['Hour', c]]
        data2 = data2.astype({c:float})
        auc = np.trapz(data2[c], x=data2.index)
        
        if plate in ABX_PLATES:
            name = columns.loc[columns.index[count/4],'name']+'_'+abx[count%4]
            type = columns.loc[columns.index[count/4],'type']
        else:
            idx = columns[columns['well']==c].index[0]
            name = columns.loc[idx,'name']
            type = columns.loc[idx,'type']
        
        
        out.append({'auc':auc,'well':c, 'name':name,'type':type})
        count+=1
    out = pd.DataFrame(out)
    
    # NOTE if there is a different control need to adjust here e.g. PM3/4, 11,12
    
    if plate in ['PM01','PM02A','PM02','PM03B']:
        negative_control = out[out['name']=='Negative Control'].index[0]
        out['auc_ratio'] = out['auc']/out.loc[negative_control,'auc'] 
    elif plate in ABX_PLATES:
        abxs = list(set([x.split('_')[0] for x in out['name'].unique().tolist()]))
        for abx in abxs:
            tmp = out[out['name'].str.contains(abx)]
            neg = tmp[tmp['name'].str.contains('1x')]
            out.loc[tmp.index,'auc_ratio'] = tmp['auc']/neg['auc'].values[0]
        
    elif plate=='PM09':
        negative_control = out[out['name']=='1% NaCl'].index[0] # may not be the best negative control
        out['auc_ratio'] = out['auc']/out.loc[negative_control,'auc'] 
    elif plate=='PM10':
        negative_control = out[out['name']=='pH 7'].index[0] # may not be the best negative control
        out['auc_ratio'] = out['auc']/out.loc[negative_control,'auc'] 
    
    if plate in ['PM04A']:
        negative_control = out[out['name']=='Negative Control'].index[0] # NOTE there are two - P source and S-source
        out['auc_ratio'] = out['auc']/out.loc[negative_control,'auc'] 
    
    return out
    
def calc_aucs(data):
    plates = data['plate_id'].unique()
    
    out = pd.DataFrame()
    for plate_id in plates:
        print('Calculating AUCs for plate_id: %s'%plate_id) # this needs to be more unique
        data2 = data[data['plate_id']==plate_id]
        type = data2['plate_type'].unique()[0]
        plate = data2[['Hour'] + WELLS] 
        
        before =  len(plate)
        plate = plate[plate['Hour']<=24]
        after = len(plate)
        print('filtering hours>24 \n\t before %i \n\t after: %i'%(before, after))
        res = calc_plate_auc(plate, type)
        res['plate_id'] = plate_id
        out = pd.concat([out,res])
        
    return out
    
def calc_fcs(data):
    plates = data['plate_id'].unique()
    
    out = pd.DataFrame()
    for plate_id in plates:
        data2 = data[data['plate_id']==plate_id]
        type = data2['plate_type'].unique()[0]
        strain = data2['strain'].unique()[0]
        plate = data2[['Hour'] + WELLS] 
        res = calc_plate_fc(plate, type)
        res['plate_id'] = plate_id
        out = pd.concat([out,res])
        
    return out
    
def growth_heatmap(out, value, columns='strain', show=True):
    res = pd.pivot_table(out, index='name', columns=columns, values=value, aggfunc=np.mean)
    
    import scipy
    sources = pd.DataFrame(scipy.stats.zscore(res), index=res.index, columns=res.columns)
    strains = pd.DataFrame(scipy.stats.zscore(res.transpose()), index=res.columns, columns=res.index)
    
    sns.clustermap(res, cmap='coolwarm')
    if show==True:
        plt.show()
    else:
        plt.savefig('binary_heatmap.svg')
        plt.savefig('binary_heatmap.png')
        plt.close()
    
    return res
    
def binary_heatmap(out, value, columns='strain', threshold=1.1, show=True):
    res = pd.pivot_table(out, index='name', columns=columns, values=value, aggfunc=np.mean)
    
    import scipy
    sources = pd.DataFrame(scipy.stats.zscore(res), index=res.index, columns=res.columns)
    strains = pd.DataFrame(scipy.stats.zscore(res.transpose()), index=res.columns, columns=res.index)
    
    res[res<threshold]=0
    res[res>=threshold]=1
    
    sums = res.transpose().sum()
    
    # only plot unique - so filter growth on all or none
    sums = sums[sums>0] # none grow
    
    sums = sums.sort_values(ascending=False)
    sums.plot(kind='bar')
    if show==True:
        plt.show()
    else:
        plt.savefig('growth_sums.svg')
        plt.savefig('growth_sum.png')
        plt.close()
        
    sums = sums[sums!=sums.max()] # all grow
    
    res_filtered = res.ix[sums.index]
    sns.clustermap(res_filtered, cmap='gray_r')
    if show==True:
        plt.show()
    else:
        plt.savefig('binary_heatmap.svg')
        plt.savefig('binary_heatmap.png')
        plt.close()
    
    return res
    
    
def PCA_plot(out, value, columns=['strain'], show=True):

    res = pd.pivot_table(out, index='name', columns=columns, values=value, aggfunc=np.mean)
    
    import scipy
    # run PCA on both strains and sources
    sources = pd.DataFrame(scipy.stats.zscore(res), index=res.index, columns=res.columns)
    strains = pd.DataFrame(scipy.stats.zscore(res.transpose()), index=res.columns, columns=res.index)
    
    names = ['strains','sources']
    idx = 0
    
    for res in [strains, sources]:
        from sklearn.decomposition import PCA
        pca = PCA(n_components=2)
        pca.fit(res.transpose().fillna(0))
        xy = pd.DataFrame(pca.components_).transpose()
        
        fig, ax = plt.subplots()
        ax.scatter(xy[0], xy[1], color='red')
        
        for i, txt in enumerate(res.index.tolist()):
            ax.annotate(txt, (xy.loc[i,0], xy.loc[i,1]))
        
        plt.xlabel('Ax1: %0.2f'%pca.explained_variance_ratio_[0])
        plt.ylabel('Ax2: %0.2f'%pca.explained_variance_ratio_[1])
        
        if show==True:
            plt.show()
        else:
            plt.savefig(names[idx]+'_PCA.svg')
            plt.savefig(names[idx]+'_PCA.png')
            plt.close()
        idx+=1
    
    
def barcharts(out, value='auc', columns=['strain'], show=True):
    res = pd.pivot_table(out, index='name', columns=columns, values=value, aggfunc=np.mean)
    
    import scipy
    sources = pd.DataFrame(scipy.stats.zscore(res), index=res.index, columns=res.columns)
    res = sources
    
    from sklearn.decomposition import PCA
    pca = PCA(n_components=2)
    pca.fit(res.transpose())
    xy = pd.DataFrame(pca.components_, columns=res.index).transpose()
    
    xy['dim1']=abs(xy[0]) 
    xy['dim2']=abs(xy[1]) 
    max_diff = xy
    for val in ['dim1','dim2']:
        max_diff = xy.sort_values(val, ascending=False)[:50]
        
        cols = 5
        rows = 10
        fig, ax_array = plt.subplots(rows, cols,squeeze=False, figsize=(12,12))
        count = 0
        for i,ax_row in enumerate(ax_array):
            for j,axes in enumerate(ax_row):
                idx = max_diff.index[count]
                tmp = out[out['name']==idx]
                mean = tmp.groupby(columns).mean() 
                std = tmp.groupby(columns).std()
                mean[value].plot(kind='bar', yerr=std, ax=axes, title='')
                axes.set_title(idx)
                if i != rows-1:
                    axes.set_xticklabels([])
                count+=1
        if show==True:
            plt.show()
        else:
            plt.savefig(val+'_barchart.svg')
            plt.savefig(val+'_barchart.png')
            plt.close()
            
def get_study_auc_matrix(session, study_name, index=['strain'], param='auc_ratio'):
    aucs = get_study_processed_data(session, study_name)
    if len(aucs)==0:
        return []
    matrix = pd.pivot_table(aucs, index=index, columns=['name','type'], values=param)
    return matrix
    
def package_study_data(session, study_data, out_folder):
    data = get_study_raw_data(session, study_data)
    for plate in data['plate_id'].unique():
        data2 = data[data['plate_id']==plate]
        del data2['media']
        data2 = data2[data2['hour']<=24].drop_duplicates()
        data2.to_csv(out_folder+'/plate_%i.csv'%plate)
    
def get_study_processed_data(session, study_name):
    metadata = session.query(Runinfo).join(Study2runinfo).join(Studies).filter(Studies.study_name==study_name)
    metadata = pd.read_sql(metadata.statement, metadata.session.bind)
    
    aucs = session.query(Growth).filter(Growth.plate_id.in_(metadata['plate_id'].tolist()))
    aucs = pd.read_sql(aucs.statement, aucs.session.bind)
    aucs_merged = pd.merge(aucs, metadata, left_on='plate_id', right_on='plate_id')
    
    plateinfo = session.query(Plateinfo).filter(Plateinfo.plate.in_(metadata['plate_type'].unique()))
    plateinfo = pd.read_sql(plateinfo.statement, plateinfo.session.bind)
    
    aucs_merged = pd.merge(aucs_merged, plateinfo, left_on=['plate_type','well'], right_on=['plate','well',])
    
    return aucs_merged
    
def get_auc_data(session, param='auc_ratio'): # get all auc data - for PCA analysis
    metadata = session.query(Runinfo).join(Study2runinfo).join(Studies)
    metadata = pd.read_sql(metadata.statement, metadata.session.bind)
    
    aucs = session.query(Growth).filter(Growth.plate_id.in_(metadata['plate_id'].tolist()))
    aucs = pd.read_sql(aucs.statement, aucs.session.bind)

    aucs = pd.merge(aucs, metadata, left_on='plate_id', right_on='plate_id')
    
    matrix = pd.pivot_table(aucs, index='strain', columns='name', values=param)
    
    return matrix
   
    

def get_study_raw_data(session, study_name, mongo_table):
    from . import io

    metadata = session.query(Runinfo).join(Study2runinfo).join(Studies).filter(Studies.study_name==study_name)
    metadata = pd.read_sql(metadata.statement, metadata.session.bind)

    data = io.get_growth_data(metadata, session, mongo_table)
    
    return data

## Adding in curve fit functions for kinetic parameters
def GompertzFunction(time,a,b,c):
    return a*np.exp(-np.exp(b-c*time))

def RichardsFunction(time,a,b,c,d,e):
    return a + (b-a)/(1+c*np.exp(-d*time))**(1/e)

def LogisticFunction(time,a,b,c):
    return a/(1+np.exp(-b*(time-c)))

def growth_correlation_map(data,plate):
    correlation_matrix = np.zeros([len(data.columns),len(data.columns)])
    for i in range(0,len(data.columns)):
        for j in range(0,len(data.columns)):
            correlation_matrix[i,j] = scipy.stats.pearsonr(data.iloc[:,i], data.iloc[:,j])[0]
    out = pd.DataFrame(correlation_matrix, columns = data.columns, index = data.columns)

    return out
