import os
basedir = os.path.abspath(os.path.dirname(__file__))


''' MONGO DB '''
MONGOPASS = '' # ASK JON FOR THIS PASSWORD

MONGOCLIENT = 'mongodb://phenommongo:%s@cluster0-shard-00-00.sauib.mongodb.net:27017,cluster0-shard-00-01.sauib.mongodb.net:27017,cluster0-shard-00-02.sauib.mongodb.net:27017/phenom?ssl=true&replicaSet=atlas-vtmenq-shard-0&authSource=admin&retryWrites=true&w=majority'%MONGOPASS

MONGO_CNX = {
    'client':MONGOCLIENT,
    'db': "phenom2",
    'collection':"growth_profiles2"
}


''' SQLLITE '''
SQLALCHEMY_DATABASE_URI = 'sqlite:///' + os.path.join(basedir, 'postil.db')
