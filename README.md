### IDMiner

"Lo que hace el programa"

##### To run this app:

You can clone or download this repo:   
```
git clone https://github.com/agustincastro/datascience-IdMiner/
```

Then cd into the repo:   
```
cd datascience-IdMiner
```
## To install dependencies using venv

Create a virtual environment called virtual_env unsing python 3 venv tool.

`apt-get install python3-venv`

`python3 -m venv vidminer`

Now create and activate a virtualenv:   
On a mac or linux:   
```
source vidminer/bin/activate
```

On a Windows:   
```
vidminer/Scripts/activate
```

Now that virtualenv is setup and active we can install the dependencies and wheels: 

Install project requirements (in some case you need to install wheels).

```
pip install wheel
pip install -r requirements.txt
```

Once the dependencies have been installed, run the application:
```
python index.py
```

Then visit http://127.0.0.1:8050/
