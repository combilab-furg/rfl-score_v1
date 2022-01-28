#-----------------------------------------------------------------------------
# Pharamcophore based SASA for complex (Author: Cheng Wang)
#-----------------------------------------------------------------------------
import os, sys, pybel
import numpy as np
import pandas as pd
#from .pharma import pharma
from .pharma import pharma

def get_columns():
    return ['S' + str(i+1) for i in range(10)]

#def runMSMS(inprot, inlig, MSMSDIR = '.'):
def runMSMS(receptor, ligand, tmp_folder_name, name):
    
    # create tmp folder for all intermediate files
    #olddir = os.getcwd()
    #os.chdir(MSMSDIR)
    #os.system('mkdir tmp')
    
    # Process Input file to be p.pdb and l.pdb
    # convert protein file to PDB if not and remov
    #ppdb = 'tmp/p.pdb'e hetatm card
    ppdb = os.path.join(tmp_folder_name, "p.pdb")
    
    __, intype = os.path.splitext(receptor)
    
    if intype[1:].lower() != 'pdb':
        #prot = pybel.readfile(intype[1:], inprot).__next__()
        prot = pybel.readfile(intype[1:], receptor).__next__()
        output = pybel.Outputfile("pdb", ppdb, overwrite=True)
        output.write(prot)
        output.close()
    else:
        # change possible HETATM to ATOM in pdb
        os.system("""sed 's/HETATM/ATOM\ \ /g' """ + receptor + " > " + ppdb)

    # convert ligand file to be PDB by openbabel
    #lpdb = 'tmp/l.pdb'
    lpdb = os.path.join(tmp_folder_name, "l.pdb")
    __, intype = os.path.splitext(ligand)
    if intype[1:].lower() != 'pdb':
        lig = pybel.readfile(intype[1:], ligand).__next__()
        output = pybel.Outputfile("pdb", lpdb, overwrite=True)
        output.write(lig)
        output.close()
    else:
        # change possible HETATM to ATOM in pdb
        os.system("""sed 's/HETATM/ATOM\ \ /g' """ + ligand + " > " + lpdb)

    #os.chdir('tmp')
    #copy atom typefiel into directory
    #os.system("cp $DXGB/atmtypenumbers .")
    # Process p.pdb/l.pdb to be p_sa.pdb/l_sa.pdb after pharma assignment
    # get full atom idx list and pharma
    ppdb2 = os.path.join(tmp_folder_name, "p_sa.pdb")
    lpdb2 = os.path.join(tmp_folder_name, "l_sa.pdb")

    pidx, ppharm = pharma(ppdb).assign(write=True, outfn = ppdb2)
    lidx, lpharm = pharma(lpdb).assign(write=True, outfn = lpdb2)

    # get subset atom idx which is nine element type
    # This have been done in pharma but still do it again
    elementint = [6, 7, 8, 9, 15, 16, 17, 35, 53]
    psub = [idx for idx in pidx if ppharm[idx][0] in elementint]
    lsub = [idx for idx in lidx if lpharm[idx][0] in elementint]

    # get element number and pharma type and assign to df1
    comp = []
    for idx in psub:
        comp.append(ppharm[idx][0:2])
    for idx in lsub:
        comp.append(lpharm[idx][0:2])

    df1 = {}
    df1['atm'] = np.array(comp)[:,0]
    df1['pharma'] = np.array(comp)[:,1]
    df1 = pd.DataFrame(df1)

    # pdb to xyzr convert
    base_path = os.path.join(__file__.rsplit("/", 1)[0], "msms")
    pdb_to_xyzr = os.path.join(base_path, "pdb_to_xyzr")
    msms = os.path.join(base_path, "msms")
    atmtypenumbers = os.path.join(base_path, "atmtypenumbers")

    p_sa_xyzr = os.path.join(tmp_folder_name, "p_sa.xyzr")
    l_sa_xyzr = os.path.join(tmp_folder_name, "l_sa.xyzr")
    pl_sa_xyzr = os.path.join(tmp_folder_name, "pl_sa.xyzr")
    
    #os.system("pdb_to_xyzr " + ppdb2 + " > p_sa.xyzr")
    command = pdb_to_xyzr + " " + ppdb2 + " " + atmtypenumbers + " > " + p_sa_xyzr
    os.system(command)

    #os.system("pdb_to_xyzr " + lpdb2 + " > l_sa.xyzr")
    command = pdb_to_xyzr + " " + lpdb2 + " " + atmtypenumbers + " > " + l_sa_xyzr
    os.system(command)
    
    #os.system("cat p_sa.xyzr l_sa.xyzr > pl_sa.xyzr")
    command = "cat " + p_sa_xyzr + " " + l_sa_xyzr + " > " +  pl_sa_xyzr
    os.system(command)

    # run msms in with radius 1.0 (if fail, will increase to be 1.1)
    p_sa_area = os.path.join(tmp_folder_name, "p_sa.area")
    l_sa_area = os.path.join(tmp_folder_name, "l_sa.area")
    pl_sa_area = os.path.join(tmp_folder_name, "pl_sa.area")
    
    log1_tmp = os.path.join(tmp_folder_name, "log1.tmp")
    log2_tmp = os.path.join(tmp_folder_name, "log2.tmp")
    log3_tmp = os.path.join(tmp_folder_name, "log3.tmp")

    #os.system("msms -if p_sa.xyzr  -af p_sa.area -probe_radius 1.0 -surface ases > log1.tmp 2>&1")
    command = msms + " -if " + p_sa_xyzr + " -af " + p_sa_area + " -probe_radius 1.0 -surface ases > " + log1_tmp + " 2>&1"    
    os.system(command)

    #os.system("msms -if l_sa.xyzr  -af l_sa.area -probe_radius 1.0 -surface ases > log2.tmp 2>&1")
    command = msms + " -if " + l_sa_xyzr + " -af " + l_sa_area + " -probe_radius 1.0 -surface ases > " + log2_tmp + " 2>&1"    
    os.system(command)

    #os.system("msms -if pl_sa.xyzr  -af pl_sa.area -probe_radius 1.0 -surface ases > log3.tmp 2>&1")
    command = msms + " -if " + pl_sa_xyzr + " -af " + pl_sa_area + " -probe_radius 1.0 -surface ases > " + log3_tmp + " 2>&1"    
    os.system(command)

    #if (os.path.isfile('p_sa.area') and os.path.isfile('l_sa.area') and os.path.isfile('pl_sa.area')) == False:
    if (os.path.isfile(p_sa_area) and os.path.isfile(l_sa_area) and os.path.isfile(pl_sa_area)) == False:
        #os.system("msms -if p_sa.xyzr  -af p_sa.area -probe_radius 1.1 -surface ases > log1.tmp 2>&1")
        command = msms + " -if " + p_sa_xyzr + " -af " + p_sa_area + " -probe_radius 1.1 -surface ases > " + log1_tmp + " 2>&1"    
        os.system(command)
    
        #os.system("msms -if l_sa.xyzr  -af l_sa.area -probe_radius 1.1 -surface ases > log2.tmp 2>&1")
        command = msms + " -if " + l_sa_xyzr + " -af " + l_sa_area + " -probe_radius 1.1 -surface ases > " + log2_tmp + " 2>&1"    
        os.system(command)

        #os.system("msms -if pl_sa.xyzr  -af pl_sa.area -probe_radius 1.1 -surface ases > log3.tmp 2>&1")
        command = msms + " -if " + pl_sa_xyzr + " -af " + pl_sa_area + " -probe_radius 1.1 -surface ases > " + log3_tmp + " 2>&1"    
        os.system(command)

        print('1.1')
    #if (os.path.isfile('p_sa.area') and os.path.isfile('l_sa.area') and os.path.isfile('pl_sa.area')) == False:
    if (os.path.isfile(p_sa_area) and os.path.isfile(l_sa_area) and os.path.isfile(pl_sa_area)) == False:
        print("SASA failed")

    # read surface area to df2
    df2 = {}
    #tmp1 = np.genfromtxt('p_sa.area', skip_header=1)[:,2]
    tmp1 = np.genfromtxt(p_sa_area, skip_header=1)[:,2]
    num_p = len(tmp1)

    #tmp2 = np.genfromtxt('l_sa.area', skip_header=1)[:,2]
    tmp2 = np.genfromtxt(l_sa_area, skip_header=1)[:,2]
    num_l = len(tmp2)
    
    #tmp3 = np.genfromtxt('pl_sa.area', skip_header=1)[:,2]
    tmp3 = np.genfromtxt(pl_sa_area, skip_header=1)[:,2]

    df2[2] = np.append(tmp1, tmp2)
    df2[3] = tmp3
    df2 = pd.DataFrame(df2)
    df = pd.concat([df1, df2], axis=1)
    df.columns = ['atm','pharma','pl','c']

    df_pro = df[0:num_p].copy()
    df_lig = df[num_p:num_p + num_l].copy()

    return df, df_pro, df_lig

def featureSASA(receptor, ligand, tmp_folder_name, name, write=False):

    # nine elements and nine pharma types
    #elemint = [6, 7, 8, 9, 15, 16, 17, 35, 53]
    #elemstr = [str(i) for i in elemint]
    pharmatype = ['P', 'N', 'DA', 'D', 'A', 'AR', 'H', 'PL', 'HA']
    outdict = {i:0 for i in pharmatype}
    outdict_pro = {i:0 for i in pharmatype}
    outdict_lig = {i:0 for i in pharmatype}

    # run MSMS
    df,df_pro,df_lig = runMSMS(receptor = receptor, ligand = ligand, tmp_folder_name = tmp_folder_name, name = name)

    ## delta SASA with clip 0 (if value less 0, cut to 0)
    df["d"] = (df["pl"] - df["c"]).clip(0,None)
    df_pro["d"] = (df_pro["pl"] - df_pro["c"]).clip(0,None)
    df_lig["d"] = (df_lig["pl"] - df_lig["c"]).clip(0,None)

    # group delta sasa by element and pharma type
    
    dfg =  df.groupby("pharma")["d"].sum()
    dfgdict =  dfg.to_dict()

    dfg_pro =  df_pro.groupby("pharma")["d"].sum()
    dfgdict_pro =  dfg_pro.to_dict()

    dfg_lig =  df_lig.groupby("pharma")["d"].sum()
    dfgdict_lig =  dfg_lig.to_dict()


    # assign grouped dict to outdict
    for i in dfgdict:
        outdict[i] = dfgdict[i]

    for i in dfgdict_pro:
        outdict_pro[i] = dfgdict_pro[i]

    for i in dfgdict_lig:
        outdict_lig[i] = dfgdict_lig[i]

    # output list
    sasalist = []
    sasalist_pro = []
    sasalist_lig = []
    for i in pharmatype:
        sasalist.append(outdict[i])
        sasalist_pro.append(outdict_pro[i])
        sasalist_lig.append(outdict_lig[i])

    sasalist.append(sum(sasalist))
    sasalist_pro.append(sum(sasalist_pro))
    sasalist_lig.append(sum(sasalist_lig))

    if write:

        sasa_dat = os.path.join(tmp_folder_name, "sasa.dat")

        f = open(sasa_dat, "w")
        f.write(" ".join([str(np.round(i,2)) for i in sasalist]) + "\n")
        f.close()

    return df, df_pro, df_lig, sasalist, sasalist_pro, sasalist_lig

def get_sasa10_descriptors(receptor, ligand, tmp_folder_name, name):
    
    rawdata, rawdata_pro, rawdata_lig, sasa, sasa_pro, sasa_lig = featureSASA(receptor = receptor, ligand = ligand, tmp_folder_name = tmp_folder_name, name = name)    
    
    # Features
    sasaTotal = sasa[-1]
    sasa_proTotal = sasa_pro[-1]
    sasa_ligTotal = sasa_lig[-1]
    sasaFeatures = sasa[0:-1]
    sasa_proFeatures = sasa_pro[0:-1]
    sasa_ligFeatures = sasa_lig[0:-1]
    
    dt = pd.DataFrame(columns = get_columns())
    dt.loc[name,:] = sasa

    dt["pdb"] = name
    dt.index = dt.pop("pdb")

    return dt
