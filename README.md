
## IDMiner


*Abstract text minnig approach to explore overrepresentated biological terms from gene lists*



### 1. Clone or download this repo

```
git clone https://github.com/sradiouy/IdMiner
```

### 2. Enter to the project directory

```
cd IdMiner
```

### 3. Create and activate a virtualenv:   


#### 3.1 On Mac or Linux:

```
apt-get install python3-venv   --- (only if needed)

python3 -m venv vidminer

source vidminer/bin/activate
````

#### 3.2 On Windows: 
 
 If it is your first time with python app, we recommend the following things: 
  
  1. installation of visual studio code (https://code.visualstudio.com/docs/python/python-tutorial)
  1. Go to https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2017 and download "Build Tools For Visual Studio" under "Tools for " Visual Studio"

Once python and visual sutdio are installed:

```
python3 -m venv vidminer

vidminer/Scripts/activate 
````
If the previous line has an error please run:

```
vidminer/Scripts/activate.bat 

```

### 4. Install project requirements (in some case you need to install wheels).

```
pip install wheel
pip install -r requirements.txt
```

### 5. Run the application.

```
python idminer.py
```

### 6. Enter to localhost in any web browser:

````
http://127.0.0.1:8050/
````


**We hope that IDMiner will be useful for your research!** :v::v::v::v:
