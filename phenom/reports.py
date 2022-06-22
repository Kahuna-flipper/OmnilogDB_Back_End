
from IPython import embed
import pandas as pd
import seaborn as sns
import pylab as plt
import numpy as np

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

def ts_plots(data, pdf, group='strain'):

    
    abx_plates = ['PM11C','PM12B']  
    
    plates = data['plate_type'].unique()    
    
    import datetime
    import numpy as np
        
    for plate in plates:
        print(plate)
        data2 = data[data['plate_type']==plate]
        data2 = data2.astype({'hour':float})
        data2 = data2.reset_index()
        
        columns = pd.read_csv('data/plate_layouts/%s.csv'%plate, header=None)
        columns = columns.rename(columns={0:'plate',1:'well',2:'name',3:'type'})
        columns.index=columns['well']
        
        max_y = 400
        
        count=0
        abx = ['1x','2x','4x','8x'] # not sure if this is correct, could be linear?
        for c in data[WELLS]:
            print(c)
            if group in METADATA:
                data3 = data2[METADATA+[c]]
            else:
                data3 = data2[METADATA+[c]+[group]]
            data3 = data3.astype({c:float})
            if plate in abx_plates:
                name = columns.loc[columns.index[count/4],'name']+'_'+abx[count%4]
                type = columns.loc[columns.index[count/4],'type']
            else:
                name = columns.loc[c,'name']
                type = columns.loc[c,'type']
            
            sns.lineplot(x='hour', y=c, data=data3, hue=group)
            plt.title(plate + ' ' + c + ' ' + name + ' (' + type +')')
            plt.ylim(0,max_y)
            # plt.show()
            
            # plt.savefig('figures/%s_%s.png'%(plate,c))
            # plt.close()
            
            pdf.savefig()  # saves the current figure into a pdf page
            plt.close()
            
            count+=1
            # '''
            


def barcharts(out, pdf, value='auc', columns=['strain']):
    res = pd.pivot_table(out, index='name', columns=columns, values=value, aggfunc=np.mean)
    
    
    for source in res.index:
        
        tmp = out[out['name']==source]
        # embed()
        mean = tmp.groupby(columns).mean() 
        std = tmp.groupby(columns).std()
        mean[value].plot(kind='bar', yerr=std, title=source)
        plt.tight_layout()
        plt.title(source)
        plt.ylabel(value)
        plt.xticks(fontsize=12, rotation=90)
        # plt.show()
        print('saving barchart for %s'%source)
        pdf.savefig()
        plt.close()
            
def PCA_plots(out, pdf, value='auc', columns=['strain'], name=''):
    res = pd.pivot_table(out, index='name', columns=columns, values=value, aggfunc=np.mean)
    import pylab as plt
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
        plt.title(name)
        pdf.savefig()
        # pdf.savefig()
        plt.close()
        idx+=1
   
    
