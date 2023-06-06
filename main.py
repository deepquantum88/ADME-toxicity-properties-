#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().system('pip install pytdc')


# In[2]:


from tdc.single_pred import ADME
from tdc.single_pred import Tox
from tdc.single_pred import HTS


# In[3]:


data=Tox(name="DILI")


# In[4]:


split=data.get_split()


# In[5]:


type(split["test"])


# In[6]:


split["test"]


# In[7]:


#split["test"].to_csv("dili_test.csv")


# In[8]:


#split["train"].to_csv("dili_train.csv")


# In[9]:


import pandas as pd
df1=pd.read_csv("dili_test.csv")
df2=pd.read_csv("dili_train.csv")


# In[10]:


data=pd.concat([df1,df2], axis=0)


# In[11]:


data


# In[12]:


from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import PandasTools
from rdkit import DataStructs
import numpy as np

def computeMorganFP(mol, depth=2, nBits=2048):
    a = np.zeros(nBits)
    try:
      DataStructs.ConvertToNumpyArray(AllChem.GetMorganFingerprintAsBitVect(mol,depth,nBits),a)
    except:
      return None
    return a

def computeRDKitFP(mol, maxLength=7, nBits=2048):
    a = np.zeros(nBits)
    try:
      DataStructs.ConvertToNumpyArray(Chem.RDKFingerprint(mol,maxPath=maxLength, fpSize=nBits),a)
    except:
      return None
    return a


# In[13]:


PandasTools.AddMoleculeColumnToFrame(frame=data, smilesCol='Drug', molCol='Molecule')
problematicSmiles = data.loc[data['Molecule'].map(lambda x: x is None)]
problematicSmiles


# In[14]:


data = data.loc[data['Molecule'].map(lambda x: x is not None)]
data.describe()


# In[15]:


data['RDKit7FP'] = data['Molecule'].map(computeRDKitFP)


# In[16]:


f=np.array([x for x in data['RDKit7FP']])


# In[17]:


y=data["Y"]


# In[18]:


len(f)


# In[19]:


f.shape


# In[ ]:




