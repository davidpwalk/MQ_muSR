# -*- coding: utf-8 -*-

# drafting template for downloading data from the msruser.psi.ch repository and evaluation

# os for checking local files/directories
import os
# pathlib for cross platform handling of paths
import pathlib
# requests to download the file
import requests
# tarfile to unpack
import tarfile
# shutil to copy the files
import shutil
# sys to exit for probably erroneous repo folder spec
import sys

import numpy as np
#import matplotlib.pyplot as plt

# outdated get_run function. this is just kept for backwards compatibility. Please use get_asy below
def get_run(run = 891, year = 2020, alpha = 2.0):
    instrument = 'GPS'
#    year = 2020
    #run = 891
    
    
    #year = 2019
    #run = 3407
    
    #%% part about getting the file from web
    storagefold = r'Python\repository'
    storagefold = pathlib.Path(storagefold)
    
    # check the local storage folder
    if not os.path.exists(storagefold):
        os.mkdir(storagefold)
        
        
    file_ID = '{:d}_{:d}_{}.txt'.format(year, run, instrument)
    
    file_fullpath = storagefold / instrument / file_ID
    
    localfile = False
    localfolder = False
    
    # check for instrument folder
    if os.path.exists(storagefold / instrument):
        localfolder = True
        if os.path.isfile(file_fullpath):
            localfile = True
    
    # download ASCII if it is not stored locally
    if not localfile:
        url = 'http://musruser.psi.ch/cgi-bin/SearchDB.cgi?AREA={}&YEAR={:d}&RUN={:d}&OUTFormat=ASCII&go=SingleFileAscii'.format(instrument, year, run)
        
    #    if requests.is_downloadable(url)
        r = requests.get(url)
#        print(r.headers.get('content-type'))
        
        # was this successful? check the content. A wrong url will give some ascii output, a proper one the Data.tgz file
        if 'Data.tgz' in r.headers.get('content-type'):
            tgzfile = 'temp.tgz' # this is a local temporary tgz file
            file = open(tgzfile,'wb')
            file.write(r.content) # which is identical to the downloaded data.tgz
            file.close()
            
            ascii_found = False
            t = tarfile.open(tgzfile)
            # it has the .ascii and the .out file. we only want the .ascii file
            mbrs = t.getmembers()
            for ii, mbr in enumerate(mbrs):
                if '.ascii' in mbr.name:
                    exfile = mbr.name
                    ascii_found = True
                    break
                
            # extract and copy it
            if ascii_found:
                t.extract(mbrs[ii])
            
                # only now create the instrument directory
                if not localfolder:
                    os.mkdir(storagefold / instrument)
                
                shutil.copy2(exfile, file_fullpath)
                
                # cleanup 
                os.remove(exfile)
            else:
                print('TODO: panic! No ascii file found. The datafile for the requested run only contained:')
                for mbr in mbrs:
                    print(mbr.name)
                    
                
            # close and remove the tarfile
            t.close()
            os.remove(tgzfile)
            
            
            
    #%% part about parsing the file
    if os.path.isfile(file_fullpath):
        # get the header first
        headerlines = 0
        header = {}
        with open(file_fullpath) as f:
            for line in f:
                if line[0] == '%':
                    headerlines+=1
                    prts = line.split(':')
                    if len(prts) == 2:
    #                    key = prts[0][1:]
                        header[prts[0][1:].strip()] = prts[1].strip()
                else:
#                    print(line)
                    break
                    
        dta_raw = np.genfromtxt(file_fullpath, skip_header=headerlines)
        
        t0_str = header['t0'].split(',')
        deltat = np.double(header['time resolution'].split(' ')[0])
        
        # offset after which data is taken
        t0_offset = 8
        
        trcs = []
        offsets = []
        trcs_x = []
        t0s = []
        for ii in range(dta_raw.shape[1]):
            
            trcs.append(dta_raw[int(t0_str[ii])+t0_offset:,ii])
            trcs_x.append(np.arange(len(trcs[ii]))*deltat)
            
            offsets.append(np.mean(dta_raw[:int(t0_str[ii])-t0_offset,ii]))
            t0s.append(int(t0_str[ii]))
            
            # keep track of smallest size
            if ii == 0:
                zero_idx = np.where(trcs[ii] == 0.0)[0]
                if len(zero_idx) > 0:
                    minlen = zero_idx[0]
                else:
                    minlen = np.min(len(trcs[ii]))
            else:
                zero_idx = np.where(trcs[ii] == 0.0)[0]
                if len(zero_idx) > 0:
                    minlen = np.min([minlen, zero_idx[0]])
                else:
                    minlen = np.min([minlen, len(trcs[ii])])
                    
        # above approach fails for datasets with too small statistics... just get everything, 
        t0s = np.array(t0s)
        offsets = np.array(offsets)
        
        minlen = dta_raw.shape[0] - np.max(t0s) - t0_offset
                
        # make more userfiendly array by cropping the endpoints
        dta = np.zeros((minlen, dta_raw.shape[1]))
        dta_x = np.arange(minlen)*deltat
        
        for ii in range(dta_raw.shape[1]):
            dta[:,ii] = trcs[ii][:minlen] - offsets[ii]
            
        # data contains the following detectors
        # 0: Forw
        # 1: Back
        # 2: Up
        # 3: Down
        # 4: Right
        # 5: Left
        
        # edit: this changed in 2023, where all is coming from root files
        dta = sort_counts(dta)
        
        
        # construct asymmetry
#        alpha = 2.0
        beta = 1.0
        asy = np.zeros((minlen, 3))
        asy[:,0] = (alpha*dta[:,1] - dta[:,0])/(alpha*beta*dta[:,1] + dta[:,0])
        asy[:,1] = (alpha*dta[:,3] - dta[:,2])/(alpha*beta*dta[:,3] + dta[:,2])
        asy[:,2] = (alpha*dta[:,5] - dta[:,4])/(alpha*beta*dta[:,5] + dta[:,4])
        
        
        return (asy, dta_x, header)
#        # some plotting
#        plt.figure(1)
#        plt.clf()
#    #    plt.plot(dta_x,dta)
#        plt.plot(dta_x,asy)

# correct implementation:
# added error and adjusted limits for background correction (further away from peak)        
def get_asy(run = 891, year = 2020, alpha = 2.0, binning = 1, start_t = 0.0, end_t = 100.0e3, single_bin = False, instrument = 'GPS'):

    # retrieve the counts
    [dta, dta_x, header] = get_counts(run, year, alpha, binning, start_t, end_t, single_bin, instrument = instrument)
            
    N_pts = dta.shape[0]
        
    # data contains the following detectors
    # 0: Forw
    # 1: Back
    # 2: Up
    # 3: Down
    # 4: Right
    # 5: Left
    
    
    # edit: this changed in 2023, where all is coming from root files
    dta = sort_counts(dta)
    
    # construct asymmetry
#        alpha = 2.0
    # alpha is differeent for each asy
    alpha_vec = [None]*3
    try:
        n_asy = len(alpha)
        
        alpha_vec[0] = alpha[0]
        alpha_vec[1] = alpha[1]
        alpha_vec[2] = alpha[2]
        
        if n_asy == 3:
            alpha_vec[2] = alpha[2]
        
    except:
        alpha_vec[0] = alpha
        alpha_vec[1] = alpha
        alpha_vec[2] = alpha
        
    beta = 1.0
    asy = np.zeros((N_pts, 3))
    asy[:,0] = (alpha_vec[0]*dta[:,1] - dta[:,0])/(alpha_vec[0]*beta*dta[:,1] + dta[:,0])
    asy[:,1] = (alpha_vec[1]*dta[:,3] - dta[:,2])/(alpha_vec[1]*beta*dta[:,3] + dta[:,2])
    # right/left detectors found to be absent for some 2014 GPS datasets
    try:
        asy[:,2] = (alpha_vec[2]*dta[:,5] - dta[:,4])/(alpha_vec[2]*beta*dta[:,5] + dta[:,4])
    except:
        pass
    
    # get error vector
    
    # that' s as implemented in the musr_TD_PSI_bin.cpp
    asy_err = np.zeros((N_pts, 3))
    asy_err[:,0] = 2*alpha_vec[0]*np.sqrt(dta[:,1]*dta[:,0]*(dta[:,1]+dta[:,0]))/(alpha_vec[0]*dta[:,1] + dta[:,0])**2
    asy_err[:,1] = 2*alpha_vec[1]*np.sqrt(dta[:,3]*dta[:,2]*(dta[:,3]+dta[:,2]))/(alpha_vec[1]*dta[:,3] + dta[:,2])**2
    try:
        asy_err[:,2] = 2*alpha_vec[2]*np.sqrt(dta[:,5]*dta[:,4]*(dta[:,5]+dta[:,4]))/(alpha_vec[1]*dta[:,5] + dta[:,4])**2
    except:
        pass
    
    # that' s from https://docs.mantidproject.org/nightly/algorithms/AsymmetryCalc-v1.html#id7
    asy_err2 = np.zeros((N_pts, 3))
    asy_err2[:,0] = np.sqrt(dta[:,0] + alpha_vec[0]**2*dta[:,1])*np.sqrt(1+asy[:,0]**2)/(dta[:,0]+alpha_vec[0]*dta[:,1])
    asy_err2[:,1] = np.sqrt(dta[:,2] + alpha_vec[1]**2*dta[:,3])*np.sqrt(1+asy[:,1]**2)/(dta[:,2]+alpha_vec[1]*dta[:,3])
    try:
        asy_err2[:,2] = np.sqrt(dta[:,4] + alpha_vec[2]**2*dta[:,5])*np.sqrt(1+asy[:,2]**2)/(dta[:,4]+alpha_vec[2]*dta[:,5])
    except:
        pass
    
    # mask out those that are zero
#        silent_ids = np.where((dta[:,1] < 0.5) | (dta[:,0] < 0.5))[0]
#        if len(silent_ids) > 0:
#            asy_err[silent_ids,0] = 1.0
#            
#        silent_ids = np.where((dta[:,3] < 0.5) | (dta[:,2] < 0.5))[0]
#        if len(silent_ids) > 0:
#            asy_err[silent_ids,1] = 1.0
#            
#        silent_ids = np.where((dta[:,5] < 0.5) | (dta[:,4] < 0.5))[0]
#        if len(silent_ids) > 0:
#            asy_err[silent_ids,2] = 1.0
    
        
    # return the asymmetry error number two. It came closer to the result of MUSRFIT than asy_err. Could be numerics, though
    return (asy, dta_x, header, asy_err2)


# function to retrieve the counts, incl eventual binning, from a dataset
def get_counts(run = 891, year = 2020, alpha = 2.0, binning = 1, start_t = 0.0, end_t = 100.0e3, single_bin = False, instrument = 'GPS'):
    
    #%% part about getting the file from web
    storagefold = r'../../Data/GPS'
    storagefold = pathlib.Path(storagefold)
    
    # check the local storage folder
    if not os.path.exists(storagefold):
        
        # hardcode here some safety... for the unlike case that someone else uses the same path, this can just be commented out...
        if (storagefold == pathlib.Path('J:/muon/pysmus/repository')): sys.exit('Seems you run this for the first time. Please set storagefold first.');
        os.mkdir(storagefold)
        
        
    file_ID = '{:d}_{:d}_{}.txt'.format(year, run, instrument)
    
    file_fullpath = storagefold / instrument / file_ID
    
    localfile = False
    localfolder = False
    
    # check for instrument folder
    if os.path.exists(storagefold / instrument):
        localfolder = True
        if os.path.isfile(file_fullpath):
            localfile = True
    
    # download ASCII if it is not stored locally
    if not localfile:
        url = 'http://musruser.psi.ch/cgi-bin/SearchDB.cgi?AREA={}&YEAR={:d}&RUN={:d}&OUTFormat=ASCII&go=SingleFileAscii'.format(instrument, year, run)
        
    #    if requests.is_downloadable(url)
        r = requests.get(url)
#        print(r.headers.get('content-type'))
        
        # was this successful? check the content. A wrong url will give some ascii output, a proper one the Data.tgz file
        if 'Data.tgz' in r.headers.get('content-type'):
            tgzfile = 'temp.tgz' # this is a local temporary tgz file
            file = open(tgzfile,'wb')
            file.write(r.content) # which is identical to the downloaded data.tgz
            file.close()
            
            ascii_found = False
            t = tarfile.open(tgzfile)
            # it has the .ascii and the .out file. we only want the .ascii file
            mbrs = t.getmembers()
            for ii, mbr in enumerate(mbrs):
                if '.ascii' in mbr.name:
                    exfile = mbr.name
                    ascii_found = True
                    break
                
            # extract and copy it
            if ascii_found:
                t.extract(mbrs[ii])
            
                # only now create the instrument directory
                if not localfolder:
                    os.mkdir(storagefold / instrument)
                
                shutil.copy2(exfile, file_fullpath)
                
                # cleanup 
                os.remove(exfile)
            else:
                print('TODO: panic! No ascii file found. The datafile for the requested run only contained:')
                for mbr in mbrs:
                    print(mbr.name)
                return []
                    
                
            # close and remove the tarfile
            t.close()
            os.remove(tgzfile)
            
            
            
    #%% part about parsing the file
    if os.path.isfile(file_fullpath):
        # get the header first
        headerlines = 0
        header = {}
        with open(file_fullpath) as f:
            for line in f:
                if line[0] == '%':
                    headerlines+=1
                    prts = line.split(':')
                    if len(prts) == 2:
    #                    key = prts[0][1:]
                        header[prts[0][1:].strip()] = prts[1].strip()
                else:
#                    print(line)
                    break
                    
        dta_raw = np.genfromtxt(file_fullpath, skip_header=headerlines)
        
        t0_str = header['t0'].split(',')
        deltat = np.double(header['time resolution'].split(' ')[0])
        
        # offset after which data is taken
        t0_offset = 8
        t0_offset = 0 # to be consistent with musrfit...
        
        # number of points before t0, where the background shall be extracted
        background_offset = 20
        
        trcs = []
        offsets = []
        trcs_x = []
        t0s = []
        for ii in range(dta_raw.shape[1]):
            
            trcs.append(dta_raw[int(t0_str[ii])+t0_offset:,ii])
            trcs_x.append(np.arange(len(trcs[ii]))*deltat)
            
            offsets.append(np.mean(dta_raw[:int(t0_str[ii])-background_offset,ii]))
            t0s.append(int(t0_str[ii]))
            
            # keep track of smallest size
            if ii == 0:
                zero_idx = np.where(trcs[ii] == 0.0)[0]
                if len(zero_idx) > 0:
                    minlen = zero_idx[0]
                else:
                    minlen = np.min(len(trcs[ii]))
            else:
                zero_idx = np.where(trcs[ii] == 0.0)[0]
                if len(zero_idx) > 0:
                    minlen = np.min([minlen, zero_idx[0]])
                else:
                    minlen = np.min([minlen, len(trcs[ii])])
                    
        # above approach fails for datasets with too small statistics... just get everything, 
        t0s = np.array(t0s)
        offsets = np.array(offsets)
        
        # print('________')
        # print(offsets)
        # print('________')
        
        minlen = dta_raw.shape[0] - np.max(t0s) - t0_offset
        
        # time axis for full dataset
        full_x = np.arange(minlen)*deltat
        
        # binning and cutting:
        cut_idx = np.where((full_x >= start_t) & (full_x <= end_t))[0]
        
        # special option to reduce a range to one single bin
        if single_bin:
            binning = cut_idx[-1] - cut_idx[0] + 1
            
            
        bin_idx = cut_idx[::binning]
            
#        if single_bin:
#            print('cutting idx from {} to {}, binning thus {}'.format(cut_idx[0],cut_idx[-1],binning) )
#            print(bin_idx)
 
        N_pts = len(bin_idx)
        
        # the cropped and binned array
        dta = np.zeros((N_pts, dta_raw.shape[1]))
        
        # center the databins
        dta_x = full_x[bin_idx] + deltat * (binning-1)/2
        
        # put data
        for ii in range(dta_raw.shape[1]):
            if binning > 0:
                for jj in range(len(bin_idx)):
                    dta[jj,ii] = np.sum(trcs[ii][bin_idx[jj]:bin_idx[jj]+binning] - offsets[ii])
            else: # directly copy if there is no binning
                dta[:,ii] = trcs[ii][bin_idx] - offsets[ii]
                
        # give back the array
        return (dta, dta_x, header)
    
    # if the file is not found, there is a problem..
    else:
        return []
    
   
# function to sort the histograms from get_counts according to the detectors of the instrument
def sort_counts(dta, run = 891, year = 2020, instrument='GPS'):
    
    
    # sorted data will the following detectors
    # 0: Forw
    # 1: Back
    # 2: Up
    # 3: Down
    # 4: Right
    # 5: Left
    dta_srt = np.zeros((dta.shape[0],6))
    
    
    # [01] Forw Run gps23_3819(11):	 8587377, 8442564
    # [02] Back Run gps23_3819(11):	 11411048, 11272124
    # [03] Up_B Run gps23_3819(11):	 6313581, 6196937
    # [04] Up_F Run gps23_3819(11):	 4031204, 3949254
    # [05] Down_B Run gps23_3819(11):	 5703204, 5599531
    # [06] Down_F Run gps23_3819(11):	 3720480, 3657936
    # [07] Right_B Run gps23_3819(11):	 1926176, 1894421
    # [08] Right_F Run gps23_3819(11):	 1361560, 1348637
    # [09] Left_B Run gps23_3819(11):	 2936547, 2890114
    # [10] Left_F Run gps23_3819(11):	 1508297, 1491197
    # [11] Mob-RL Run gps23_3819(11):	 2194020, 2154797
    # TODO: we are ignoring Mob-RL!!!
    
    
    """
    # that is the 'old' PSIbin detector arrangement
    for ii in range(6):
        dta_srt[:,ii] = dta[:,ii]
    """
    
    if (dta.shape[1] == 11):
        print('Grouping detectors with FWDVETO')
        # fwd/bwd
        dta_srt[:,0] = dta[:,0]
        dta_srt[:,1] = dta[:,1]
        
        # up
        dta_srt[:,2] = dta[:,2] + dta[:,3]
        # down
        dta_srt[:,3] = dta[:,4] + dta[:,5]
        # right
        dta_srt[:,4] = dta[:,6] + dta[:,7]
        # left
        dta_srt[:,5] = dta[:,8] + dta[:,9]
    
    
    elif (dta.shape[1] == 12):
        print('Grouping detectors without FWDVETO')
        # fwd/bwd
        dta_srt[:,0] = dta[:,0] + dta[:,1]
        dta_srt[:,1] = dta[:,2]
        
        # up
        dta_srt[:,2] = dta[:,3] + dta[:,4]
        # down
        dta_srt[:,3] = dta[:,5] + dta[:,6]
        # right
        dta_srt[:,4] = dta[:,7] + dta[:,8]
        # left
        dta_srt[:,5] = dta[:,9] + dta[:,10]
        
    elif (dta.shape[1] == 6):
        print('1:1 detector mapping (bin files)')
        for ii in range(6):
            dta_srt[:,ii] = dta[:,ii]
    elif (dta.shape[1] == 4):
        print('1:1 detector mapping with four detectors')
        for ii in range(4):
            dta_srt[:,ii] = dta[:,ii]
    elif (dta.shape[1] == 2):
        print('1:1 detector mapping with two detectors')
        for ii in range(2):
            dta_srt[:,ii] = dta[:,ii]
            # SUPERHACK: here is the flame arrangement, just based on the number of detectors found... ugly, but good for a weekends beamtime...
    elif (dta.shape[1] == 8):
        print('FLAME detector mapping')
        # fwd/bwd
        dta_srt[:,0] = dta[:,0]
        dta_srt[:,1] = dta[:,1]
        
        # left
        dta_srt[:,2] = dta[:,3] + dta[:,6] + dta[:,7]
        
        # right
        dta_srt[:,3] = dta[:,2] + dta[:,4] + dta[:,5]
        
    else:
        print('Found {:} histograms in ascii file, do not know what to do with that...'.format(dta.shape[1]))
    
    
    
    # run and year info kept for the moment, but probably not needed, unless one goes beyond GPS...
    
    return dta_srt
    
    
   
    
# experimental difference function
# used to extract the difference of the decay asymmetry. For instance uwave On/Off, or phase 0/180
def get_asy_diff(runs = [891, 892], year = 2020, alpha = 2.0, binning = 1, start_t = 0.0, end_t = 100.0e3, single_bin = False):

    if len(runs) != 2:
        print('difference function requires two runs. You specified {} runs'.format(len(runs)))
        
    # retrieve the counts
    [d1, d1_x, header1] = get_counts(runs[0], year, alpha, binning, start_t, end_t, single_bin)
    N_pts1 = d1.shape[0]
    
    # edit: this changed in 2023, where all is coming from root files
    d1 = sort_counts(d1)
    
    [d2, d2_x, header2] = get_counts(runs[1], year, alpha, binning, start_t, end_t, single_bin)
    N_pts2 = d2.shape[0]
    
    # edit: this changed in 2023, where all is coming from root files
    d2 = sort_counts(d2)
    
    # combine the datasets together
    if N_pts1 != N_pts2:
        print('The number of points of the two datasets is different ({} vs {}). This should not happen with the same binnig parameters'.format(N_pts1,N_pts2))
    
    dta_x = d1_x.copy()
    
    # build the difference
    dta_nom = d1 - d2
    dta_denom = d1 + d2
    
    N_pts = N_pts1
    
    # this will be used for the error bar. Double check if that is correct
    dta = dta_denom
    
                
    # data contains the following detectors
    # 0: Forw
    # 1: Back
    # 2: Up
    # 3: Down
    # 4: Right
    # 5: Left
    
    # construct asymmetry
#        alpha = 2.0
    # alpha is differeent for each asy
    alpha_vec = [None]*3
    try:
        n_asy = len(alpha)
        
        alpha_vec[0] = alpha[0]
        alpha_vec[1] = alpha[1]
        alpha_vec[2] = alpha[2]
        
        if n_asy == 3:
            alpha_vec[2] = alpha[2]
        
    except:
        alpha_vec[0] = alpha
        alpha_vec[1] = alpha
        alpha_vec[2] = alpha
        
    beta = 1.0
    asy = np.zeros((N_pts, 3))
    asy[:,0] = (alpha_vec[0]*dta_nom[:,1] - dta_nom[:,0])/(alpha_vec[0]*beta*dta_denom[:,1] + dta_denom[:,0])
    asy[:,1] = (alpha_vec[1]*dta_nom[:,3] - dta_nom[:,2])/(alpha_vec[1]*beta*dta_denom[:,3] + dta_denom[:,2])
    asy[:,2] = (alpha_vec[2]*dta_nom[:,5] - dta_nom[:,4])/(alpha_vec[2]*beta*dta_denom[:,5] + dta_denom[:,4])
    
    # get error vector
    
    # that' s as implemented in the musr_TD_PSI_bin.cpp
    asy_err = np.zeros((N_pts, 3))
    asy_err[:,0] = 2*alpha_vec[0]*np.sqrt(dta[:,1]*dta[:,0]*(dta[:,1]+dta[:,0]))/(alpha_vec[0]*dta[:,1] + dta[:,0])**2
    asy_err[:,1] = 2*alpha_vec[1]*np.sqrt(dta[:,3]*dta[:,2]*(dta[:,3]+dta[:,2]))/(alpha_vec[1]*dta[:,3] + dta[:,2])**2
    asy_err[:,2] = 2*alpha_vec[2]*np.sqrt(dta[:,5]*dta[:,4]*(dta[:,5]+dta[:,4]))/(alpha_vec[1]*dta[:,5] + dta[:,4])**2
    
    
    # that' s from https://docs.mantidproject.org/nightly/algorithms/AsymmetryCalc-v1.html#id7
    asy_err2 = np.zeros((N_pts, 3))
    asy_err2[:,0] = np.sqrt(dta[:,0] + alpha_vec[0]**2*dta[:,1])*np.sqrt(1+asy[:,0]**2)/(dta[:,0]+alpha_vec[0]*dta[:,1])
    asy_err2[:,1] = np.sqrt(dta[:,2] + alpha_vec[1]**2*dta[:,3])*np.sqrt(1+asy[:,1]**2)/(dta[:,2]+alpha_vec[1]*dta[:,3])
    asy_err2[:,2] = np.sqrt(dta[:,4] + alpha_vec[2]**2*dta[:,5])*np.sqrt(1+asy[:,2]**2)/(dta[:,4]+alpha_vec[2]*dta[:,5])
    
    # mask out those that are zero
#        silent_ids = np.where((dta[:,1] < 0.5) | (dta[:,0] < 0.5))[0]
#        if len(silent_ids) > 0:
#            asy_err[silent_ids,0] = 1.0
#            
#        silent_ids = np.where((dta[:,3] < 0.5) | (dta[:,2] < 0.5))[0]
#        if len(silent_ids) > 0:
#            asy_err[silent_ids,1] = 1.0
#            
#        silent_ids = np.where((dta[:,5] < 0.5) | (dta[:,4] < 0.5))[0]
#        if len(silent_ids) > 0:
#            asy_err[silent_ids,2] = 1.0
    
        
    # return the asymmetry error number two. It came closer to the result of MUSRFIT than asy_err. Could be numerics, though
    return (asy, dta_x, header1, asy_err2)