### IDMiner

"Lo que hace el programa"

##### To run IDMiner:

You can __clone__ or __download__ this repo:   
```
git clone https://github.com/sradiouy/IdMiner
```

Then cd into the repo:   
```
cd IdMiner
```
## To install dependencies using venv

1. Now create and activate a virtualenv:   


### 1.1 On Mac or Linux:

```
apt-get install python3-venv   --- (only if needed)

python3 -m venv vidminer

source vidminer/bin/activate
````

### 1.2 On Windows: 
 
 If it is your first time with python app, we recommend the following things: 
  
  1. installation of visual studio code (https://code.visualstudio.com/docs/python/python-tutorial)
  2. Go to https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2017 and download "Build Tools For Visual Studio" under "Tools for " Visual Studio"

Once python and visual sutdio are installed:

```
python3 -m venv vidminer

vidminer/Scripts/activate 
````
If the previous line has an error please run:

```
vidminer/Scripts/activate 

```

Now that virtualenv is setup and active we can install the dependencies and wheels: 

### 2. Install project requirements (in some case you need to install wheels).

```
pip install wheel
pip install -r requirements.txt
```

Once the dependencies have been installed, run the application:
```
python idminer.py
```

Then visit http://127.0.0.1:8050/


**We hope that IDMiner will be useful for your research!** :v::v::v::v:
