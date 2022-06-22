# BiologDB

This repo serves as the backend for BiologDB

It loads curated biolog data in csv format into the BiologDB

It also performs simple growth curve analyses like AUC and growth fold changes

Basic reporting functionality is also included



# Getting started

1. Copy and rename the config.example.py file to config.py

2. Update the MONGO_CNX line with information from Jon (best practice not to save login info in repos)

3. There are two options for installing needed dependencies:
    1. Create a new conda env with a fresh python install and all needed packages: ```conda env create -f biologdb.yml```
    2. Use your exsiting python install and pip to install requirements: ```pip install -r requirements.txt```

4. If using conda, activate the new env with ```conda activate biologdb```

5. Test the three sample scripts:

    1. example_load_runs.py

    2. example_queries.py

    3. example_analysis.py