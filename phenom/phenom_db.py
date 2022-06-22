import sqlalchemy
from IPython import embed

from sqlalchemy.ext.declarative import declarative_base

from sqlalchemy import create_engine, MetaData, TEXT, Integer, Table, Column, ForeignKey, Float, Boolean

from sqlalchemy.orm import relationship

Base = declarative_base()

class Runinfo(Base):
    
    __tablename__ = 'runinfo'
    
    plate_id = Column(Integer, primary_key=True, autoincrement=True)
    study = Column(TEXT, nullable=False) # MAY MOVE THIS TO SEPARATE TABLE
    plate_type= Column(TEXT, nullable=False)
    species = Column(TEXT, nullable=False)
    strain = Column(TEXT, nullable=False)
    media = Column(TEXT, nullable=False)
    rep = Column(TEXT, nullable=False)
    setup_time = Column(TEXT, nullable=False)
    
    # clusters = relationship("GeneCluster", back_populates="strain")
    # special_genes = relationship("SpGene", back_populates="strain")
    # genes = relationship("Gene", back_populates="strain")
    
    def __repr__(self):
        return "<Plate(plate_id='%s', plate_type='%s', strain='%s', media='%s'>" % (self.plate_id, self.plate_type, self.strain, self.media)
        
class Studies(Base):
    
    __tablename__ = 'studies'
    
    study_id = Column(Integer, primary_key=True, autoincrement=True)
    study_name = Column(TEXT, nullable=False)
    study_description= Column(TEXT, nullable=True)
    pc1 = Column(Float, nullable=True)
    pc2 = Column(Float, nullable=True)
    pc3 = Column(Float, nullable=True)
    comparator = Column(TEXT, nullable=True)
    
    def __repr__(self):
        return "<Study(study_id='%s', study_name='%s'>" % (self.study_id, self.study_name)
        
class Study2runinfo(Base):

    __tablename__ = 'study2runinfo'
    
    s2r_id = Column(Integer, primary_key=True, autoincrement=True)
    study_id = Column(Integer, ForeignKey('studies.study_id'),nullable=False)
    plate_id = Column(Integer, ForeignKey('runinfo.plate_id'), nullable=False)
    
    
# '''
class Plateinfo(Base):

    __tablename__ = 'plateinfo'

    id = Column(Integer, primary_key=True, autoincrement=True)
    
    # TO match all growth results on plate & well
    
    plate= Column(TEXT, nullable=False)
    well= Column(TEXT, nullable=False)
    cas = Column(TEXT, nullable=True)
    kegg= Column(TEXT, nullable=True)
    name= Column(TEXT, nullable=False)
    type= Column(TEXT, nullable=False)
    
    def __repr__(self):
        return "<PlateInfo(id='%s'" % (self.id)
        
class Pca(Base):

    __tablename__ = 'pca'
    
    id = Column(Integer, primary_key=True, autoincrement=True)
    
    # TO match all growth results on plate & well
    dim1 = Column(Float, nullable=False)
    dim2 = Column(Float, nullable=False)
    dim3 = Column(Float, nullable=False)
    type = Column(TEXT, nullable=False)
    name = Column(TEXT, nullable=False)
    comparator = Column(TEXT, nullable=False)
    study_id = Column(Integer, ForeignKey('studies.study_id'),nullable=False)
    
    def __repr__(self):
        return "<PCA(id='%s'" % (self.id)
        
        
# '''
        
class Growth(Base):
    
    __tablename__ = 'growth'

    auc_id = Column(Integer, primary_key=True, autoincrement=True)
    auc = Column(Float, nullable=False)
    # name= Column(TEXT, nullable=False)
    # type= Column(TEXT, nullable=False)
    well= Column(TEXT, nullable=False)
    # kegg= Column(TEXT, nullable=True)
    # cas= Column(TEXT, nullable=True)
    auc_ratio= Column(Float, nullable=False)
    plate_id = Column(Integer, ForeignKey('runinfo.plate_id'), nullable=False)
    fold_change = Column(Float, nullable=True)
    log2FC = Column(Float, nullable=True)
    
    def __repr__(self):
        return "<AUC(plate_id='%s', auc='%s'>" % (self.plate_id, self.auc)
        
'''        
class Growth(Base):
    
    __tablename__ = 'growth'
    
    id = Column(Integer, primary_key=True, autoincrement=True)
    plate_id = Column(Integer, ForeignKey('runinfo.plate_id'), nullable=False)
    hour = Column(Float, nullable=False)
    A01 = Column(Float, nullable=False)
    A02 = Column(Float, nullable=False)
    A03 = Column(Float, nullable=False)
    A04 = Column(Float, nullable=False)
    A05 = Column(Float, nullable=False)
    A06 = Column(Float, nullable=False)
    A07 = Column(Float, nullable=False)
    A08 = Column(Float, nullable=False)
    A09 = Column(Float, nullable=False)
    A10 = Column(Float, nullable=False)
    A11 = Column(Float, nullable=False)
    A12 = Column(Float, nullable=False)
    B01 = Column(Float, nullable=False)
    B02 = Column(Float, nullable=False)
    B03 = Column(Float, nullable=False)
    B04 = Column(Float, nullable=False)
    B05 = Column(Float, nullable=False)
    B06 = Column(Float, nullable=False)
    B07 = Column(Float, nullable=False)
    B08 = Column(Float, nullable=False)
    B09 = Column(Float, nullable=False)
    B10 = Column(Float, nullable=False)
    B11 = Column(Float, nullable=False)
    B12 = Column(Float, nullable=False)
    C01 = Column(Float, nullable=False)
    C02 = Column(Float, nullable=False)
    C03 = Column(Float, nullable=False)
    C04 = Column(Float, nullable=False)
    C05 = Column(Float, nullable=False)
    C06 = Column(Float, nullable=False)
    C07 = Column(Float, nullable=False)
    C08 = Column(Float, nullable=False)
    C09 = Column(Float, nullable=False)
    C10 = Column(Float, nullable=False)
    C11 = Column(Float, nullable=False)
    C12 = Column(Float, nullable=False)
    D01 = Column(Float, nullable=False)
    D02 = Column(Float, nullable=False)
    D03 = Column(Float, nullable=False)
    D04 = Column(Float, nullable=False)
    D05 = Column(Float, nullable=False)
    D06 = Column(Float, nullable=False)
    D07 = Column(Float, nullable=False)
    D08 = Column(Float, nullable=False)
    D09 = Column(Float, nullable=False)
    D10 = Column(Float, nullable=False)
    D11 = Column(Float, nullable=False)
    D12 = Column(Float, nullable=False)
    E01 = Column(Float, nullable=False)
    E02 = Column(Float, nullable=False)
    E03 = Column(Float, nullable=False)
    E04 = Column(Float, nullable=False)
    E05 = Column(Float, nullable=False)
    E06 = Column(Float, nullable=False)
    E07 = Column(Float, nullable=False)
    E08 = Column(Float, nullable=False)
    E09 = Column(Float, nullable=False)
    E10 = Column(Float, nullable=False)
    E11 = Column(Float, nullable=False)
    E12 = Column(Float, nullable=False)
    F01 = Column(Float, nullable=False)
    F02 = Column(Float, nullable=False)
    F03 = Column(Float, nullable=False)
    F04 = Column(Float, nullable=False)
    F05 = Column(Float, nullable=False)
    F06 = Column(Float, nullable=False)
    F07 = Column(Float, nullable=False)
    F08 = Column(Float, nullable=False)
    F09 = Column(Float, nullable=False)
    F10 = Column(Float, nullable=False)
    F11 = Column(Float, nullable=False)
    F12 = Column(Float, nullable=False)
    G01 = Column(Float, nullable=False)
    G02 = Column(Float, nullable=False)
    G03 = Column(Float, nullable=False)
    G04 = Column(Float, nullable=False)
    G05 = Column(Float, nullable=False)
    G06 = Column(Float, nullable=False)
    G07 = Column(Float, nullable=False)
    G08 = Column(Float, nullable=False)
    G09 = Column(Float, nullable=False)
    G10 = Column(Float, nullable=False)
    G11 = Column(Float, nullable=False)
    G12 = Column(Float, nullable=False)
    H01 = Column(Float, nullable=False)
    H02 = Column(Float, nullable=False)
    H03 = Column(Float, nullable=False)
    H04 = Column(Float, nullable=False)
    H05 = Column(Float, nullable=False)
    H06 = Column(Float, nullable=False)
    H07 = Column(Float, nullable=False)
    H08 = Column(Float, nullable=False)
    H09 = Column(Float, nullable=False)
    H10 = Column(Float, nullable=False)
    H11 = Column(Float, nullable=False)
    H12 = Column(Float, nullable=False)
    
    # run = relationship("GeneCluster", back_populates="strain")
    # special_genes = relationship("SpGene", back_populates="strain")
    # genes = relationship("Gene", back_populates="strain")
    
    def __repr__(self):
        return "<Measurement(id='%s', plate_id='%s', hour='%s'>" % (self.id, self.plate_id, self.hour)
'''

