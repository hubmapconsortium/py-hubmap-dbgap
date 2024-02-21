#!/usr/bin/env python
# coding: utf-8

# # Daily Report
# This script can be turned into an AirFlow DAG.

# In[ ]:


import hubmapbags

token = '<this-is-my-token>'

hubmapbags.utilities.clean()
df = hubmapbags.reports.daily(token=token)

