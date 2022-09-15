#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import csv


# In[61]:


list1=np.loadtxt("parsnp.vcf.csv", skiprows=1,delimiter=',', dtype=str)
mut=list1[:,0:1]
nt=list1[:,1:2]
ref=list1[:,2:3]
lfv=list1[:,3:4]


# In[62]:


list2=np.loadtxt("primer_regions.csv", skiprows=1,delimiter=',', dtype=str)
primername=list2[:,0:1]
fasta=list2[:,1:2]
sense=list2[:,2:3]
start=list2[:,3:4]
stop=list2[:,4:5]


# In[63]:


report=[]
for x in range(len(primername)):
    for y in range(len(list1)):
        p=primername[x]
        s=int(start[x])
        t=int(stop[x])
        r=ref[y]
        l=lfv[y]
        m=int(mut[y])
        if t>=m>=s:
            a=(p,r,l,m,s,t)
            report.append(a)


# In[69]:


list1=np.loadtxt("Baylor+mason.csv", skiprows=1,delimiter=',', dtype=str)
mut=list1[:,0:1]
ref=list1[:,1:2]
lfv=list1[:,2:3]


# In[71]:


report=[]
for x in range(len(primername)):
    for y in range(len(list1)):
        p=primername[x]
        s=int(start[x])
        t=int(stop[x])
        r=ref[y]
        l=lfv[y]
        m=int(mut[y])
        if t>=m>=s:
            a=(p,r,l,m,s,t)
            report.append(a)


# In[68]:


len(report)

