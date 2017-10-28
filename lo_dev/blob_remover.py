#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code to delete selected blobs or containers
"""

import os
import sys
pth = os.path.abspath('../../LiveOcean/alpha')
if pth not in sys.path:
    sys.path.append(pth)
import Lfun
Ldir = Lfun.Lstart()

from datetime import datetime, timedelta
from azure.storage.blob import BlockBlobService
from azure.storage.blob import PublicAccess

# ****************** CASE-SPECIFIC CODE *****************

from datetime import datetime
start_time = datetime.now()

dt0 = datetime(2017,1,1) # start time
dt1 = datetime(2017,7,31) # end time

dt = dt0
while dt <= dt1:
    Ldir['date_string'] = dt.strftime('%Y.%m.%d')
    print('Deleting Azure files for ' + Ldir['date_string'])
    f_string = 'f' + Ldir['date_string']
    # Azure commands
    ff_string = f_string.replace('.','') # azure does not like dots in container names
    # account name and key
    azu_dict = Lfun.csv_to_dict(Ldir['data'] + 'accounts/azure_pm_2015.05.25.csv')
    account = azu_dict['account']
    key = azu_dict['key']
    containername = ff_string
    # get a handle to the account
    blob_service = BlockBlobService(account_name=account, account_key=key)
    blob_service.create_container(containername)
    blob_service.set_container_acl(containername, public_access=PublicAccess.Container)
    result = blob_service.delete_container(containername)
    print(' - result = ' + str(result))
    dt = dt + timedelta(days=1)
print('DONE')